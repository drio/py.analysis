package main

import (
	"fmt"
	"github.com/drio/drio.go/common/files"
	"strings"
	"io"
	"bufio"
	"os"
	"strconv"
	"net/http"
	"io/ioutil"
	"encoding/json"
)

type genome map[string][]int16

func (g genome) add(chrm string, coor int) {
	if _, present := g[chrm]; !present {
		g[chrm] = make([]int16, 1000000)
	}
	// We may have to grow the slice
	if coor >= len(g[chrm]) {
		t := make([]int16, coor + 10000000)
		for i := range g[chrm] {
    	t[i] = g[chrm][i]
		}
		g[chrm] = t
	}
	g[chrm][coor]++
}

func (g genome) get(chrm string, coor int) int16 {
	if _, present := g[chrm]; !present {
		return 0
	}
	if coor >= len(g[chrm]) {
		return 0
	}
	return g[chrm][coor]
}

type bedReader struct {
	fd *os.File
	reader *bufio.Reader
	coor, end int // current coordinate, end for current bed entry
	chrm string // current chrm
}

func (br *bedReader) enumSites() (string, int, bool) {
	var err error
	var start, end string

	if br.coor == 0 || br.coor > br.end {
		var line []byte
		line, err = br.reader.ReadSlice('\n')
		if err != nil {
			if err == io.EOF {
				return "", 0, true
			} else {
				panic(err)
			}
		}
		sLine := strings.Split(string(line), "\t")
		br.chrm, start, end = sLine[0], sLine[1], strings.Trim(sLine[2], "\n")
		if br.coor, err = strconv.Atoi(start); err != nil {
			panic(err)
		}
		if br.end, err = strconv.Atoi(end); err != nil {
			panic(err)
		}
	}

	br.coor++;
  return br.chrm, br.coor-1, false
}

type query struct {
	Sites []site
}

type site struct {
  Chrm  string
  Start int
}

type answer []int16

func genHandler(gnm genome) http.HandlerFunc {
	return func (w http.ResponseWriter, r *http.Request) {
		var body []byte
		var err error

		if body, err = ioutil.ReadAll(r.Body); err != nil {
			panic(err)
		}
		q := &query{}
		if err = json.Unmarshal(body, q); err != nil {
			panic(err)
		}

		//fmt.Println("Query: " , q)
		aAnswer := make(answer, len(q.Sites))
		for i:=0; i<len(aAnswer); i++ {
			aAnswer[i] = gnm.get(q.Sites[i].Chrm, q.Sites[i].Start)
		}
		wire, _ := json.Marshal(aAnswer)
		fmt.Fprintf(w, string(wire))
	}
}

func startServer(gnm genome) {
  fmt.Println("Starting server ...")
  http.HandleFunc("/", genHandler(gnm))
	if err := http.ListenAndServe(":8080", nil); err != nil {
		fmt.Println(err)
	}
}

func main() {
	fd, reader := files.Xopen("-");
	defer fd.Close()
	br := &bedReader{fd, reader, 0, 0, ""}
	g  := make(genome)
	count := 0
	for chrm, coor, done := br.enumSites(); !done; chrm, coor, done = br.enumSites() {
		count++
		g.add(chrm, coor)
		if count == 100000 {
			fmt.Println(chrm, coor)
			count = 0
		}
	}
	fmt.Println("Done.")
	startServer(g)
}
