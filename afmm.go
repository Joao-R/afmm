package afmm

import (
	"container/heap"
	"container/list"
	"image"
	"math"
	"sync"
)

type DataGrid struct {
	U      []float64
	T      []float64
	f      []uint8
	colNum int
	rowNum int
}

type pixelHeap struct {
	data      *DataGrid
	heapIndex map[int]int
	pixelIds  []int
}

func (h pixelHeap) Len() int { return len(h.pixelIds) }
func (h pixelHeap) Less(i, j int) bool {
	return h.data.T[h.pixelIds[i]] < h.data.T[h.pixelIds[j]]
}

func (h pixelHeap) Swap(i, j int) {
	h.pixelIds[i], h.pixelIds[j] = h.pixelIds[j], h.pixelIds[i]
	h.heapIndex[h.pixelIds[i]] = i
	h.heapIndex[h.pixelIds[j]] = j
}

func (h *pixelHeap) Push(x any) {
	if h.data.f[x.(int)] == 1 {
		h.heapIndex[x.(int)] = h.Len()
		h.pixelIds = append(h.pixelIds, x.(int))
		return
	}

	heap.Fix(h, h.heapIndex[x.(int)])
}

func (h *pixelHeap) Pop() any {
	x := h.pixelIds[h.Len()-1]
	h.pixelIds = h.pixelIds[:h.Len()-1]
	return x
}

func vonNeumannNeighborhood(idx int, colNum int) [4]int {
	x := idx % colNum
	y := idx / colNum
	y = y * colNum

	return [4]int{y + x - 1, y - colNum + x, y + x + 1, y + colNum + x}
}

func mooreNeighborhood(idx int, colNum int) [8]int {
	x := idx % colNum
	y := idx / colNum
	y = y * colNum
	ym1 := y - colNum
	yp1 := y + colNum

	return [8]int{ // CCW around pixel x,y starting from x-1,y-1
		ym1 + x - 1,
		ym1 + x,
		ym1 + x + 1,
		y + x + 1,
		yp1 + x + 1,
		yp1 + x,
		yp1 + x - 1,
		y + x - 1,
	}
}

// We use this to ensure the initial pixel in boundary
// detection lies outside the object being skeletonized
func (d *DataGrid) safeMooreNeighborhood(idx int) [8]int {
	neighbors := mooreNeighborhood(idx, d.colNum)
	var offset int
	for i, neighbor := range neighbors {
		if d.f[neighbor] == 0 {
			offset = i
			break
		}
	}
	var neighborsStartingOutside [8]int

	for i := range neighbors {
		neighborsStartingOutside[i] = neighbors[(i+offset)%8]
	}

	return neighborsStartingOutside
}

func (imgData *DataGrid) ParseImage(img *image.Image) {
	bounds := (*img).Bounds()
	imgData.colNum = bounds.Max.X - bounds.Min.X + 2
	imgData.rowNum = bounds.Max.Y - bounds.Min.Y + 2

	imgData.f = make([]uint8, imgData.colNum*imgData.rowNum)
	imgData.T = make([]float64, imgData.colNum*imgData.rowNum)

	var xx, yy int // relative x and y
	yy = 1
	for y := bounds.Min.Y; y < bounds.Max.Y; y++ {
		xx = 1
		for x := bounds.Min.X; x < bounds.Max.X; x++ {

			r, g, b, _ := (*img).At(x, y).RGBA()
			lum := 0.299*float64(r) + 0.587*float64(g) + 0.114*float64(b)

			if lum > 32768 {
				imgData.f[xx+yy*imgData.colNum] = 1
				imgData.T[xx+yy*imgData.colNum] = math.MaxFloat64
			} else {
				imgData.f[xx+yy*imgData.colNum] = 0
				imgData.T[xx+yy*imgData.colNum] = 0
			}

			xx++
		}
		yy++
	}
}

func (state *DataGrid) initFMM(band *pixelHeap) {
	var idx int
	/* Boundary detection */
	band.data = state
	band.heapIndex = make(map[int]int, (state.rowNum)*(state.colNum)/2) // not used for FMM only, but needed since Pop and Push interact with it
	for y := 1; y < state.rowNum-1; y++ {
		for x := 1; x < state.colNum-1; x++ {
			idx = y*state.colNum + x
			if state.f[idx] == 1 {
				for _, j := range vonNeumannNeighborhood(idx, state.colNum) {
					if state.f[j] == 0 {
						band.heapIndex[idx] = len(band.pixelIds)
						band.pixelIds = append(band.pixelIds, idx)
						state.T[idx] = 0
						state.f[idx] = 2 // band
						break
					}
				}
			}
		}
	}

	heap.Init(band)
}

func (state *DataGrid) initAFMM(band *pixelHeap, startInFront bool) {
	var idx int
	/* Boundary detection */

	idxToList := make(map[int]*list.Element)
	bandList := list.New()

	addToList := bandList.PushBack
	if !startInFront {
		addToList = bandList.PushFront
	}

	band.data = state
	band.heapIndex = make(map[int]int, (state.rowNum)*(state.colNum)/2)
	for y := 1; y < state.rowNum-1; y++ {
		for x := 1; x < state.colNum-1; x++ {
			idx = y*state.colNum + x
			if state.f[idx] == 1 {
				for _, j := range vonNeumannNeighborhood(idx, state.colNum) {
					if state.f[j] == 0 {

						idxToList[idx] = addToList(idx)
						band.pixelIds = append(band.pixelIds, idx)
						state.T[idx] = 0
						state.f[idx] = 3 // 3: band uninitialized
						break
					}
				}
			}
		}
	}

	/* Initialize U */

	var found bool
	var current int
	count := 0

	for bandList.Len() > 0 {

		current = bandList.Remove(bandList.Front()).(int)
		state.U[current] = float64(count)
		state.f[current] = 2
		count++

		/* propagation */
		found = true
		for found {
			found = false
			neighbors := state.safeMooreNeighborhood(current)
			for _, j := range neighbors {
				if state.f[j] == 3 {
					current = j
					state.U[current] = float64(count)
					state.f[current] = 2
					count++

					bandList.Remove(idxToList[current])

					found = true
					break
				}
			}
		}
	}

	heap.Init(band)
}

func solve(idx1 int, idx2 int, T []float64, f []uint8, solution *float64) {
	var r, s float64

	if f[idx1] == 0 {
		if f[idx2] == 0 {
			r = math.Sqrt(2 - ((T[idx1] - T[idx2]) * (T[idx1] - T[idx2])))
			s = (T[idx1] + T[idx2] - r) * .5
			if s >= T[idx1] && s >= T[idx2] {
				if s < *solution {
					*solution = s
				} else {
					s += r
					if s >= T[idx1] && s >= T[idx2] {
						if s < *solution {
							*solution = s
						}
					}
				}
			}
		} else {
			if 1+T[idx1] < *solution {
				*solution = 1 + T[idx1]
			}
		}
	} else {
		if f[idx2] == 0 {
			if 1+T[idx2] < *solution {
				*solution = 1 + T[idx2]
			}
		}
	}
}

func propagateU(idx int, U []float64, f []uint8, neighbors [4]int) {
	var a, m, M float64
	var counter float64
	a = 0
	m = math.MaxFloat64
	M = -1.0
	counter = 0.0
	for _, idx := range neighbors {
		if f[idx] == 0 {
			a += U[idx]
			if U[idx] < m {
				m = U[idx]
			}
			if U[idx] > M {
				M = U[idx]
			}
			counter++
		}
	}
	a /= counter

	if (M - m) < 2.0 {
		U[idx] = a
	}
}

func (d *DataGrid) stepAFMM(band *pixelHeap) {
	var solution float64

	current := heap.Pop(band).(int)

	d.f[current] = 0
	for _, neighbor := range vonNeumannNeighborhood(current, d.colNum) {
		if d.f[neighbor] != 0 {

			otherNeighbors := vonNeumannNeighborhood(neighbor, d.colNum)

			solution = (*d).T[neighbor]

			solve(otherNeighbors[0], otherNeighbors[1], d.T, d.f, &solution)
			solve(otherNeighbors[2], otherNeighbors[1], d.T, d.f, &solution)
			solve(otherNeighbors[0], otherNeighbors[3], d.T, d.f, &solution)
			solve(otherNeighbors[2], otherNeighbors[3], d.T, d.f, &solution)

			d.T[neighbor] = solution

			heap.Push(band, neighbor)

			if d.f[neighbor] == 1 {
				d.f[neighbor] = 2
				d.U[neighbor] = d.U[current]
				propagateU(neighbor, d.U, d.f, otherNeighbors) // mooreNeighborhood(neighbor, d.colNum, true))
			}
		}
	}
}

func (d *DataGrid) stepFMM(band *pixelHeap) {
	var solution float64

	current := heap.Pop(band).(int)

	d.f[current] = 0
	for _, neighbor := range vonNeumannNeighborhood(current, d.colNum) {
		if d.f[neighbor] != 0 {

			otherNeighbors := vonNeumannNeighborhood(neighbor, d.colNum)

			solution = d.T[neighbor]

			solve(otherNeighbors[0], otherNeighbors[1], d.T, d.f, &solution)
			solve(otherNeighbors[2], otherNeighbors[1], d.T, d.f, &solution)
			solve(otherNeighbors[0], otherNeighbors[3], d.T, d.f, &solution)
			solve(otherNeighbors[2], otherNeighbors[3], d.T, d.f, &solution)

			d.T[neighbor] = solution

			heap.Push(band, neighbor)

			if d.f[neighbor] == 1 {
				d.f[neighbor] = 2
			}

		}
	}
}

func (state *DataGrid) afmm(startInFront bool) {
	var band pixelHeap
	state.initAFMM(&band, startInFront)

	for band.Len() > 0 {
		state.stepAFMM(&band)
	}
}

func FMM(img *image.Image) []float64 {
	var state DataGrid
	state.ParseImage(img)

	var band pixelHeap
	state.initFMM(&band)

	for band.Len() > 0 {
		state.stepFMM(&band)
	}

	var newIdx, oldIdx int
	DT := make([]float64, (state.colNum-2)*(state.rowNum-2))
	for y := 1; y < state.rowNum-1; y++ {
		for x := 1; x < state.colNum-1; x++ {
			oldIdx = y*state.colNum + x
			newIdx = (y-1)*(state.colNum-2) + (x - 1)
			DT[newIdx] = state.T[oldIdx]
		}
	}

	return DT
}

func AFMM(img *image.Image) []float64 {
	var stateFirst DataGrid
	stateFirst.ParseImage(img)

	stateFirst.U = make([]float64, stateFirst.colNum*stateFirst.rowNum)

	mask := make([]uint8, stateFirst.colNum*stateFirst.rowNum)
	copy(mask, stateFirst.f)

	var stateLast DataGrid
	stateLast.colNum, stateLast.rowNum = stateFirst.colNum, stateFirst.rowNum
	stateLast.f = make([]uint8, stateFirst.colNum*stateFirst.rowNum)
	stateLast.T = make([]float64, stateFirst.colNum*stateFirst.rowNum)
	stateLast.U = make([]float64, stateFirst.colNum*stateFirst.rowNum)
	copy(stateLast.f, stateFirst.f)
	copy(stateLast.T, stateFirst.T)
	copy(stateLast.U, stateFirst.U)

	var wg sync.WaitGroup

	wg.Add(1)
	go func() {
		defer wg.Done()
		stateFirst.afmm(true)
	}()

	wg.Add(1)
	go func() {
		defer wg.Done()
		stateLast.afmm(false)
	}()

	wg.Wait()

	gradU := make([]float64, (stateFirst.colNum-2)*(stateFirst.rowNum-2))

	var newIdx, oldIdx int
	var deltaUFirst, deltaULast float64
	for y := 1; y < stateFirst.rowNum-1; y++ {
		for x := 1; x < stateFirst.colNum-1; x++ {
			oldIdx = y*stateFirst.colNum + x
			if mask[oldIdx] == 0 {
				continue
			}
			newIdx = (y-1)*(stateFirst.colNum-2) + (x - 1)

			deltaUFirst = 0
			deltaULast = 0

			for _, neighbor := range mooreNeighborhood(oldIdx, stateFirst.colNum) {

				if mask[neighbor] == 0 {
					continue
				}

				difference := math.Abs(stateFirst.U[neighbor] - stateFirst.U[oldIdx])
				if deltaUFirst < difference {
					deltaUFirst = difference
				}

				difference = math.Abs(stateLast.U[neighbor] - stateLast.U[oldIdx])
				if deltaULast < difference {
					deltaULast = difference
				}
			}

			gradU[newIdx] = deltaUFirst
			if deltaULast < deltaUFirst {
				gradU[newIdx] = deltaULast
			}

			if gradU[newIdx] <= 3 {
				gradU[newIdx] = 0
			}
		}
	}

	return gradU
}
