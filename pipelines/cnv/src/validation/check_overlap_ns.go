// vim: set ts=2 noet sw=2:

/*
 * Compares two fasta files (one chrm each) and reports
 * bp locations where N is present in both files:
 */
package main

import (
	"bufio"
	"fmt"
	"github.com/drio/drio.go/common/files"
	"os"
	"runtime"
)

type Hits struct {
	bf    []bool // "bit" field
	first bool
	Count int
}

func (h *Hits) setup() {
	h.bf = make([]bool, 250000000)
	h.first = true
	h.Count = 0
}

func (h *Hits) update(bf *bufio.Reader, chrm string) int {
	mem := runtime.MemStats{}
	p := 0
	process_chrm := false
	num_ns := 0

	done := func() {
		fmt.Fprintf(os.Stderr, "\n")
		h.first = false
		return num_ns
	}

	defer done()

	for line := range files.IterLines(bf) {
		if line[0] == '>' {
			if process_chrm { // we have already process our chrm, get out
				return -1
			}
			if line[1:] == chrm {
				process_chrm = true
			}
			continue
		}

		if process_chrm {
			for _, r := range line {
				c := string(r)
				if c != "\n" {
					if c == "N" {
						num_ns++
						if h.first {
							h.bf[p] = true
						} else {
							if h.bf[p] {
								h.Count += 1
							}
						}
					}
					p++
					if p%1000000 == 0 {
						runtime.ReadMemStats(&mem)
						fmt.Fprintf(os.Stderr, "%d; [TotalAlloc: %f Mb]; Ns:%d\r", p, float64(mem.TotalAlloc/1000000), num_ns)
					}
				}
			}
		}
	}
}

func main() {
	if len(os.Args) != 5 {
		fmt.Println("Usage: tool <file1.fa> <file2.fa> <chrm_first> <chrm_second>")
		os.Exit(1)
	}

	hits := Hits{}
	hits.setup()

	chrm := os.Args[3]
	f1, bf1 := files.Xopen(os.Args[1])
	defer f1.Close()
	fmt.Fprintf(os.Stderr, "Loading 1st file... (chrm: %s)\n", chrm)
	ns_in_first := hits.update(bf1, chrm)

	chrm = os.Args[4]
	f2, bf2 := files.Xopen(os.Args[2])
	defer f2.Close()
	fmt.Fprintf(os.Stderr, "Loading 2st file... (chrm: %s)\n", chrm)
	ns_in_second := hits.update(bf2, chrm)

	fmt.Fprintf(os.Stdout, "ns_in_first %s ns_in_second %s overlap %s", ns_in_first, ns_in_second, hits.Count)
}
