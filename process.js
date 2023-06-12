const fs = require('fs');
const os = require('os');
const path = require('path');

const ROMPTH = /^OMP_NUM_THREADS=(\d+)/;
const RGRAPH = /^Loading graph .*\/(.*?)\.mtx \.\.\./m;
const RORDER = /^order: (\d+) size: (\d+) (?:\[\w+\] )?\{\}/m;
const RRESLT = /^\{-(.+?)\/\+(.+?) \[(\d+)\] batch, (.+?) threads\} -> \{(.+?)\/(.+?)ms, (.+?) iters, (.+?) passes, (.+?) modularity\} (.+)/m;




// *-FILE
// ------

function readFile(pth) {
  var d = fs.readFileSync(pth, 'utf8');
  return d.replace(/\r?\n/g, '\n');
}

function writeFile(pth, d) {
  d = d.replace(/\r?\n/g, os.EOL);
  fs.writeFileSync(pth, d);
}




// *-CSV
// -----

function writeCsv(pth, rows) {
  var cols = Object.keys(rows[0]);
  var a = cols.join()+'\n';
  for (var r of rows)
    a += [...Object.values(r)].map(v => `"${v}"`).join()+'\n';
  writeFile(pth, a);
}




// *-LOG
// -----

function readLogLine(ln, data, state) {
  ln = ln.replace(/^\d+-\d+-\d+ \d+:\d+:\d+ /, '');
  if (ROMPTH.test(ln)) {
    var [, omp_num_threads] = ROMPTH.exec(ln);
    state.omp_num_threads   = parseFloat(omp_num_threads);
  }
  else if (RGRAPH.test(ln)) {
    var [, graph] = RGRAPH.exec(ln);
    if (!data.has(graph)) data.set(graph, []);
    state.graph = graph;
  }
  else if (RORDER.test(ln)) {
    var [, order, size] = RORDER.exec(ln);
    state.order = parseFloat(order);
    state.size  = parseFloat(size);
  }
  else if (RRESLT.test(ln)) {
    var [, batch_deletions_size, batch_insertions_size, batch_index, num_threads, preprocessing_time, time, iterations, passes, modularity, technique] = RRESLT.exec(ln);
    data.get(state.graph).push(Object.assign({}, state, {
      batch_deletions_size:  parseFloat(batch_deletions_size),
      batch_insertions_size: parseFloat(batch_insertions_size),
      batch_index: parseFloat(batch_index),
      num_threads: parseFloat(num_threads),
      preprocessing_time:    parseFloat(preprocessing_time),
      time:        parseFloat(time),
      iterations:  parseFloat(iterations),
      passes:      parseFloat(passes),
      modularity:  parseFloat(modularity),
      technique,
    }));
  }
  return state;
}

function readLog(pth) {
  var text  = readFile(pth);
  var lines = text.split('\n');
  var data  = new Map();
  var state = {};
  for (var ln of lines)
    state = readLogLine(ln, data, state);
  return data;
}




// PROCESS-*
// ---------

function processCsv(data) {
  var a = [];
  for (var rows of data.values()) {
    for (var row of rows)
      a.push(row);
  }
  return a;
}




// MAIN
// ----

function main(cmd, log, out) {
  var data = readLog(log);
  if (path.extname(out)==='') cmd += '-dir';
  switch (cmd) {
    case 'csv':
      var rows = processCsv(data);
      writeCsv(out, rows);
      break;
    case 'csv-dir':
      for (var [graph, rows] of data)
        writeCsv(path.join(out, graph+'.csv'), rows);
      break;
    default:
      console.error(`error: "${cmd}"?`);
      break;
  }
}
main(...process.argv.slice(2));
