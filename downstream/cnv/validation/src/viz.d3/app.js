var app = function() {
  var csv = 'data/all.csv',
      data,
      dotplot;

  data = function(rs) {
    var _data = {},
        rows;

    function load() {
      rows = rs;
      return _data;
    }

    _data.to_Array = function(min, prop) {
      var a = [];
      for (var i=0; i<rows.length; ++i)
        if (+rows[i].min_size == min)
          a.push([rows[i].threshold, rows[i][prop]]);
      return a;
    };

    return load();
  };

  dotplot = function(data, div_id, h, w, xlabel, ylabel, title) {
    var dp = {},
        svg, g,
        r = 3, padding_x = 45, padding_y = 20,
        max_x = d3.max(data, function(d) { return d[0]; }),
        max_y = d3.max(data, function(d) { return d[1]; }),

        x = d3.scale.linear().domain([0, max_x]).range([0, w]),
        y = d3.scale.linear().domain([max_y, 0]).range([0, h]),

        x_axis = d3.svg.axis().scale(x).orient("bottom"),
        y_axis = d3.svg.axis().scale(y).orient("left");

    function create() {
      svg = d3.select(div_id)
                  .append("svg")
                  .attr('heigth', h)
                  .attr('weigth', w);

      svg.append("text")
          .attr("class", "label")
          .attr("x", (w - (w/2)))
          .attr("y", 0)
          .attr("dy", ".75em")
          .text(title)
          .attr("class", "title");

      svg.append("text")
          .attr("class", "label")
          .attr("text-anchor", "end")
          .attr("y", 0)
          .attr("x", -(h - (h/2)))
          .attr("dy", ".75em")
          .attr("transform", "rotate(-90)")
          .text(ylabel);

      svg.append("text")
          .attr("class", "label")
          .attr("x", (w - (w/2)))
          .attr("y", h + padding_x)
          .attr("dy", ".75em")
          .text(xlabel);

      svg.append("g")
         .attr("class", "axis")
         .attr("transform", "translate(" + padding_x + "," + (h + padding_y) + ")")
         .call(x_axis);

      svg.append("g")
         .attr("class", "axis")
         .attr("transform", "translate(" + padding_x + "," + padding_y + ")")
         .call(y_axis);

      g = svg.append("g")
         .attr("transform", "translate(" + padding_x + "," + padding_y + ")");

      g.selectAll("circle")
          .data(data)
          .enter()
          .append("circle")
          .attr("cx", function (d) { return x(d[0]); })
          .attr("cy", function (d) { return y(d[1]); })
          .attr("r", r);

      return dp;
    }

    return create();
  };


  d3.csv(csv, function(rows) {
    var x = 0, y = 1, min_size = 1000,
        d = data(rows),
        dp, sens;

    [ "one", "two" ].forEach(function(v) {
      d3.select("body").append("div").attr("id", v);
    });

    dotplot(d.to_Array(min_size, "sensitivity"), "#one", 150, 400, "threshold", "sensitivity", "min event size = 1000");
    dotplot(d.to_Array(min_size, "fdr"), "#two", 150, 400, "threshold", "FDR", "min event size = 1000");

  });
}();
