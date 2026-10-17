# chemrust-nasl

Bond-length-constrained new available single-atom sites searching engine for any given structures
Written in `rust` for speed and stability. Part of the `chemrust` toolchain.

## Introduction

This program is designed to carry out a thorough search of all possible bonding location for a new **single** atom to the atoms of any 2D or 3D material models with **a given target bondlength**, and to generate the calculation seed files accordingly. This program aims to conduct a rational, rapid, consistent and comprehensive model building workflows to free the computational chemistry scientists from the laborious efforts of creating possible adsorbate-surface interaction models by hands.
The program therefore naturally fits the demand of big data analysis in computational chemistry, as those massive amounts of data source from the massive amounts of material model inputs!

## Highlights

- Both the goal and the algorithm is original and unique.
- Mathematically powered to generate the positions with high accuracy and consistency, away from the painful efforts of building numerous models by hand.
- Easily find the multi-coordinated positions under the given bondlength,
- Generate the possible results of singly/doubly-coordinated sites without unwanted interactions with atoms around.
- It is especially suitable for researching the mounting sites of adsorbates on a large surface of a material, since symmetry is not considered and the spatial effect from the surface size is not affected.
- Designed to welcome less experienced computer users as well as the "veterans"."
- Can be used online in a browser or offline in terminal.
- Easily adapted to shell scripts to execute batch tasks.

## Example

```terminal
$ chemrust-nasl-app
> Filepath of the model cell file: ./scanner_test_models/black_phosphorus.cell
> Element symbol of the new atom:  78 Pt
> What is the target bondlength (Å)? 2.9
> Enter the fractional coordinate range to search in the direction of x-axis: 0
> Enter the fractional coordinate range to search in the direction of x-axis: 1
> Enter the fractional coordinate range to search in the direction of y-axis: 0
> Enter the fractional coordinate range to search in the direction of y-axis: 1
> Enter the fractional coordinate range to search in the direction of z-axis: 0
> Enter the fractional coordinate range to search in the direction of z-axis: 1
> Please name the directory for exported seed files:  Pt_2.9_black_phosphorus
> Please specify the location of castep psuedopotentials directory:  e.g.: /home/user/Potentials
> Quality for k-point sampling? Coarse
> Use edft method for the `metals_method` option in CASTEP?(y/n or yes/no) No
> Run mode of program Fast

Special multi-coordinated sites search completed.
Found 186 multi-coordinated positions;
Found 157 possible doubly-coordinated positions;
Found 100 possible singly-coordinated positions;
Results have been written to Pt_2.9_black_phosphorus
```

## Installation

Download the compiled binary from our website.

- supported platform:
  - x86_64-unknown-linux-gnu (dynamically linked, hence the host machine need to install the required dependencies if necessary, e.g. `glibc`. Further readings are in the [docs]().)
  - x86_64-unknown-linux-musl
  - macOS
  - Windows

## Usage

The program is currently offered as a CLI (command-line interface) program, which required to be run in a terminal emulator with any proper shell.
If it is installed to locations that are included in your shell's environmental variable list, you can directly call `chemrust-nasl-app` in the command-line. Or, specify the path to it, `/PATH/TO/THE/BINARY/chemrust-nasl-app`.
The program currently offers two running modes: an interactive mode and a read-from-config mode. First of all you can invoke the program by `chemrust-nasl-app -h` to get the basic help information:

```terminal
New Adsorption Site Locator (NASL) for catalyst materials. Written in rust. Part of the `chemrust` toolchain. This is the binary crate.

Usage: chemrust-nasl-app [OPTIONS] [CONFIG_LOC]

Arguments:
  [CONFIG_LOC]

Options:
  -m, --mode <MODE>  [possible values: c, i]
  -h, --help         Print help (see more with '--help')
  -V, --version      Print version
```

The interactive mode will guide you step by step, from giving the input model file (supports `.cell` from `castep` and `.cif` of the International Union of Crystallography), specifying the desired new element and bondlength, to the remaining necessary settings to generate the seed files for computation job submissions.

You will need to provide (also see the above **Example** section):

1. Path to the model. E.g. `./graphene.cif`, `./graphdiyne.cell`.
2. Element symbol of the new atom. E.g. `Cu`.
3. Target bondlength in Å. E.g. 2.2
4. The search range along the `x,y,z` axes by the fractional coordinates, from 0.0 to 1.0. Default to `0.0 - 1.0` for searching in the whole lattice.
5. Path to the exported files. If the path does not exist the program will create it for you automatically.
6. Specify the quality of k-point sampling. Default is `Coarse`.
7. Use `edft` or `dm` as the electronic minimization method in castep.
   It is recommended that use `edft` for rare-earth elements involved models, while `dm` is enough for models without rare-earth elements..
8. Running modes of the program. `Fast` means generating the resulted models files and castep job files without the relatively time-consuming process of copying the psuedopotentials files to the destination. `Full` would copy the psuedopotentials. `Post` is you can execute the copying after `Fast` has been done. The others are for debug use.

You can enter the interactive mode by `chemrust-nasl-app` or `chemrust-nasl-app -m i`.

The read-from-config mode (`-m c` or `--mode c`) is meant for faster execution and/or batch processing. It will read from a `yaml` format file which contains all the required items with the predesignated format offered by us.

An example `config.yaml`:

```yaml
model_path: ./scanner_test_models/from_mingzi/NiFeLDH_tem1.cell
new_element: Cu
target_bondlength: 2.2
x_range: [0.0, 1.0]
y_range: [0.0, 1.0]
z_range: [0.0, 1.0]
export_dir: demo/NiFeLDH_tem1_Cu_2.2
kpoint_quality: Coarse
edft: false
```

Suppose you have the required file `config.yaml` in the current directory. Run the program in read-from-config mode as follows:

```
chemrust-nasl-app -m c config.yaml
```

The program will start immediately, if the `config.yaml` is properly written without issues.
