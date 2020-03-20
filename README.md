
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HESmanip

The goal of HESmanip is to …

## Example

This is a basic example which shows you how to solve a common problem:

## Docker

This package has been developed in docker based on the
`rocker/tidyverse` image, to access the development environment enter
the following at the command line (with an active docker daemon
running),

``` bash
docker pull n8thangreen/hesmanip:0.1
docker run --rm -p 8787:8787 -e USER=HESmanip -e PASSWORD=HESmanip n8thangreen/hesmanip:0.1
```

The rstudio client can be accessed on port `8787` at `localhost` (or
your machines ip). The default username is HESmanip and the default
password is HESmanip. It can be found on Docker Hub here
<https://hub.docker.com/repository/docker/n8thangreen/hesmanip>.
