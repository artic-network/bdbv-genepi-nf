# Build stage: create the conda environment and strip what is not needed at runtime
FROM condaforge/miniforge3:latest AS build

COPY environment.yml /tmp/environment.yml

RUN mamba env create -f /tmp/environment.yml -p /opt/env \
    && mamba clean -afy \
    && find /opt/env -type f -name '*.a' -delete

# Runtime stage: slim base with only the environment copied over (no conda, no package cache)
FROM debian:bookworm-slim

# procps provides `ps`, which Nextflow uses to collect task metrics
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/env /opt/env
COPY bin/ /opt/bdbv-pipeline/bin/

ENV PATH=/opt/env/bin:/opt/bdbv-pipeline/bin:$PATH

ENV MPLCONFIGDIR="."

CMD ["/bin/bash"]
