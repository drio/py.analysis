// vim: set ts=2 noet sw=2:

/*
 * Compares two fasta files (one chrm each) and reports
 * bp locations where N is present in both files:
*
 * Example:
 * 234 2
 * 111 1
 * ...
 * In bp position 234 both files have a N
 * In bp position 111 only one of the files has the N
 *
 * Quick way to check the number of Ns in a file:
 * $ cat chr21.fa.kmer.masked | ruby -ane 'BEGIN{@c=0}; next if $_[0] == ">"; $_.each_char {|c| @c +=1 if c == "N"}; END{puts @c}'
 * 29852808
 * If we do
 * $ go run $SRC/validation/check_overlap_ns.go chr21.fa.kmer.masked chr21.fa.kmer.masked 50000000
 * we should get 29852808
 */
package main

import (
	"os"
	"fmt"
	"github.com/drio/drio.go/common/files"
	"bufio"
	"runtime"
)

type Hits struct {
	bf []bool // "bit" field
	first bool
	Count int
}

func (h *Hits) setup() {
	h.bf = make([]bool, 250000000)
	h.first = true
	h.Count = 0
}

func (h *Hits) update(bf *bufio.Reader, chrm string)  {
	mem := runtime.MemStats{}
	p := 0
	process_chrm := false

	done := func() {
		fmt.Fprintf(os.Stderr, "\n")
		h.first = false
	}

	defer done()

	for line:= range files.IterLines(bf) {
  	if line[0] == '>' {
			if process_chrm { // we have already process our chrm, get out
				return
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
						if h.first {
							h.bf[p] = true
						} else {
							if h.bf[p] {
								h.Count += 1
							}
						}
					}
					p++
					if (p % 1000000 == 0) {
						runtime.ReadMemStats(&mem)
						fmt.Fprintf(os.Stderr, "%d; [TotalAlloc: %f Mb]\r", p, float64(mem.TotalAlloc/1000000))
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
	hits.update(bf1, chrm)

	chrm = os.Args[4]
	f2, bf2 := files.Xopen(os.Args[2])
	defer f2.Close()
	fmt.Fprintf(os.Stderr, "Loading 2st file... (chrm: %s)\n", chrm)
	hits.update(bf2, chrm)

	fmt.Println(hits.Count)
}
