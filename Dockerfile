##TORUN:
# cd "C:\Users\Nathan\Documents\R\HESmanip"
# docker container ls  
# docker kill 1c87f9c37e6a 
# docker build -t my-r-image .
# docker run -e USER=HESmanip -e PASSWORD=HESmanip --rm -p 8787:8787 my-r-image

FROM rocker/verse:3.3.2

MAINTAINER "Nathan Green" ngreen1@ic.ac.uk

RUN mkdir /home/rstudio/HESmanip

#RUN R -e "install.packages('gapminder', repos = 'http://cran.us.r-project.org')"

ADD . /home/rstudio/HESmanip

RUN R -e "devtools::install('home/rstudio/HESmanip')"