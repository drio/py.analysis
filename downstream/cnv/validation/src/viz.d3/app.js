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

    _data.mins = function() {
      var a = {};
      for (var i=0; i<rows.length; ++i)
        a[rows[i].min_size] = true;
      return d3.keys(a);
    };

    return load();
  };

  dotplot = function(data, div_id, h, w, xlabel, ylabel, title_text) {
    var dp = {},
        svg, g, title,
        r = 4, padding_x = 45, padding_y = 20,
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

      title = svg.append("text")
          .attr("x", (w - (w/2)))
          .attr("y", 0)
          .attr("dy", ".75em")
          .text(title_text)
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

    dp.update = function(new_data, new_min) {
      g.selectAll("circle")
          .data(new_data)
          .transition()
          .duration(3000)
          .delay(0)
          .attr("cx", function (d) { return x(d[0]); })
          .attr("cy", function (d) { return y(d[1]); })
          .attr("r", r);

      title.text(new_min);
    };

    return create();
  };


  // Main
  d3.csv(csv, function(rows) {
    var d = data(rows), m = 1000,
        h = 200, w = 400,
        dp_one, dp_two;

    var select = d3.select("#text").append("select");

    select.selectAll("option")
      .data(d.mins())
      .enter()
      .append("option")
      .text(function(d) { return d;})
      .attr("value", function(d) {return d;});

    d3.selectAll("select").on("change", function() {
      var new_min = this.value;
      dp_one.update(d.to_Array(new_min, "sensitivity"), new_min);
      dp_two.update(d.to_Array(new_min, "fdr"), new_min);
    });

    dp_one = dotplot(d.to_Array(m, "sensitivity"), "#two", h, w, "threshold", "sensitivity", m);
    dp_two = dotplot(d.to_Array(m, "fdr"), "#two", h, w, "threshold", "FDR", m);
  });

}();
