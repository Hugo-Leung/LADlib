# LAD analysis containers

This repository contains Apptainer definition files and build scripts for an
AlmaLinux 9 container with HCANA and LADlib. The optional GEANT4 definition
adds CLHEP and GEANT4 on top of the base image.

## Requirements

- [Apptainer](https://apptainer.org/) to build from a `.def` file or to build
  the Dockerfile as an Apptainer image.
- [Docker](https://www.docker.com/) or [Podman](https://podman.io/) to build
  and run an OCI image from the Dockerfile.
- Network access during builds to download packages and source code

## Build the base image

### Using `al9.def`

From the repository directory, build the Apptainer image:

```sh
apptainer build lad-al9.sif al9.def
```

The definition installs the packages listed in `packages`, builds HCANA from
the `develop` branch.

### Using `Dockerfile`

The base image can also be built from `Dockerfile` with Docker or Podman:

```sh
docker build -t lad-al9 .
# or
podman build -t lad-al9 .
```

Or use the Dockerfile as an Apptainer build source:

```sh
apptainer build lad-al9.sif dockerfile:.
```

The Dockerfile follows the same general install flow as `al9.def`, but uses
Docker instructions in place of Apptainer sections: `ENV` sets the runtime
environment, `ENTRYPOINT` invokes `hcana` (like `%runscript`), and build-time
`RUN` commands create the directories and perform the HCANA/LADlib checks.
The checks run during the Docker build, serving a similar purpose to `%test`.

The Dockerfile also defines `ARG CACHEBUST` before the package-install and
build step. Pass it a new value to invalidate that layer and rebuild the
dependencies and software:

```sh
docker build --build-arg CACHEBUST="$(date +%s)" -t lad-al9 .
# or
podman build --build-arg CACHEBUST="$(date +%s)" -t lad-al9 .
# or, when building the Dockerfile with Apptainer
apptainer build --build-arg CACHEBUST="$(date +%s)" lad-al9.sif dockerfile:.
```

The Dockerfile currently sets only the `DEF_HASH` image label. Unlike the
Apptainer definition, it does not set `LADLIB_HASH`, `HCANA_HASH`, or
`ROOT_VERSION`: Dockerfile `LABEL` values cannot evaluate shell commands such
as `git rev-parse` or `root-config --version` during the build. Those values
can instead be recorded in a file in the image if needed.

## Add GEANT4

`geant4.def` uses `lad-al9.sif` as its base image. With that image in the
current directory, build the extended image:

```sh
apptainer build lad-al9-geant4.sif geant4.def
```

This layer builds CLHEP 2.4.7.1 and GEANT4 v11.4.0, including GEANT4 data.

## Use an image

Start an interactive shell:

```sh
apptainer shell lad-al9.sif
```

On ifarm, bind `/w` and `/cache` to make those host directories available
inside the container:

```sh
apptainer shell -B /w,/cache lad-al9.sif
```

Run a command in the container:

```sh
apptainer exec -B /w,/cache lad-al9.sif <command>
```

In the base image, `%runscript` invokes `hcana` and forwards any arguments
passed to `apptainer run`:

```sh
apptainer run -B /w,/cache lad-al9.sif [arguments...] <replay script>
```

## Repository files

- `al9.def` and the shell scripts build the base HCANA/LADlib image.
- `geant4.def`, `clhep.sh`, and `geant4.sh` build the optional GEANT4 layer.
- `packages` lists the AlmaLinux packages installed in the base image.
