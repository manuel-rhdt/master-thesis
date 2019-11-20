# Gillespie

## Installation

You have to install the dependencies from `requirements.txt` first.

## Usage

Create `configuration.toml` to configure the individual actions. An example
configuration is provided in `configuration_default.toml`.

Then from the project directory just run

```sh
python3 mutual_information
```

from the shell to execute the calculation.

The configuration file format used is [TOML](https://github.com/toml-lang/toml) whichs
is designed to be both human-readable and easy to process.

### Advanced configuration

The configuration file can be used as a [Jinja2](https://palletsprojects.com/p/jinja/) Template.
