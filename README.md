# SemAgentsToolkit
## Setup

```
(v1.6) pkg> activate .
(SemAgentsToolkit) pkg> instantiate
```

## Usage
Two example simulations are demonstrated in `scripts`

## Evluation
First, generate the Pluto notebook:
```
generateNotebook("any_name.jl", "/path/to/the/csv/files")
```

Then, open that notebook:
```
using Pluto
Pluto.run(notebook="any_name.jl")
```

