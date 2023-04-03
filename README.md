# ACS-1-armageddon

<div align=center>
<img src="https://github.com/Gekkoli/The-Hazard-of-small-Asteroids--package/blob/main/img-folder/astronomerss.jpg/">
</div>

## Installation

To install the module and any pre-requisites, from the base directory run
```
conda env create -f environment.yml
```  

## Downloading postcode data

The post code data can be downloaded by executing the python scipt. 
To download the postcode data
```
python download_data.py
```

## Automated testing

To run the pytest test suite, from the base directory run
```
pytest tests/
```

## Documentation

To generate the documentation (in html format)
```
python -m sphinx docs html
```

See the `docs` directory for the generated documentation

## Example usage

For example usage see `example.py` in the examples folder:
```
python examples/example.py
```

## More information

For more information on the developed tool, see `AirburstSolverUsage.ipynb` and `DamageMapperUsage.ipynb`.
