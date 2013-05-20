package main

import (
	"fmt"
	"github.com/drio/drio.go/common/files"
	"strings"
	"os"
	"strconv"
)

const UNDER = "_"
const DEF_MIN_NUM_SAMPLES_WITH_SNP = 1

func main() {

	minNumSamples := DEF_MIN_NUM_SAMPLES_WITH_SNP
	if len(os.Args) > 1 {
		if i, err := strconv.ParseInt(os.Args[1], 10, 32); err != nil {
			panic(err)
		} else {
			minNumSamples = int(i)
		}
	}

	fp, r := files.Xopen("-")
	defer fp.Close()

	s_ids := []string{}
	counts := make(map[string]int)

	prev_pos := ""
	for l := range(files.IterLines(r)) {
		s   := strings.Split(l, " ")
		id  := s[0]
		pos := s[1] + UNDER + s[2]

		if prev_pos != pos {
			if prev_pos != "" && len(s_ids) > 1 {
				for _, id := range(s_ids) {
					counts[id]+=1
				}
			}
			prev_pos = pos
			s_ids = s_ids[:0]
		}

		s_ids = append(s_ids, id)
	}

	if len(s_ids) > minNumSamples {
		for _, id := range(s_ids) {
			counts[id]+=1
		}
	}

	gen_size := float32(33986384)
	for id, c := range(counts) {
		fmt.Println(id, ":", (float32(c)/gen_size)*1000)
	}
}
