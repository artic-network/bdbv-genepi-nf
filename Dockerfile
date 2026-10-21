FROM condaforge/mambaforge:latest AS conda

COPY . ./bdbv-pipeline/

RUN /opt/conda/bin/mamba env create -f ./bdbv-pipeline/environment.yml

ENV PATH=/opt/conda/envs/bdbv-pipeline/bin:$PATH

ENV MPLCONFIGDIR="."

CMD ["/bin/bash"]
