# Single Benchmark Execution

## List of Benchmarks

```@eval
using EarthBox
import Markdown
options = BenchmarksManager.get_options()
liststr = join(["- `$(options[id].option_name)`" for id in keys(options)], "\n")
Markdown.parse(liststr)
```
!!! tip "quick search for benchmark description"
    Highlight a benchmark name in the list above and use `Ctl-F` + `Enter`
    to navigate to a detailed description in the **Benchmark Descriptions** section below.

## Benchmark Runner

```@docs
BenchmarksManager.run_benchmark
```