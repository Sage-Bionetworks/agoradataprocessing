# ⚠️ DEPRECATED REPOSITORY

**This repository has been deprecated and is no longer maintained.**
- No further development will occur in this repository.
- This repository will not be monitored.


## New Repository Location

Agora's data processing is now managed via the agora-data-tools project:

**[https://github.com/Sage-Bionetworks/agora-data-tools](https://github.com/Sage-Bionetworks/agora-data-tools)**


For information, see the [agora-data-tools README](https://github.com/Sage-Bionetworks/agora-data-tools/blob/dev/README.md)).

----

[![docker-cloud-automated](https://img.shields.io/docker/cloud/automated/sagebionetworks/agoradataprocessing.svg)](https://cloud.docker.com/u/sagebionetworks/repository/docker/sagebionetworks/agoradataprocessing) [![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/sagebionetworks/agoradataprocessing.svg)](https://cloud.docker.com/u/sagebionetworks/repository/docker/sagebionetworks/agoradataprocessing/builds)
# Agora Data Processing

R package and command line utility for processing data from Synapse to Agora. This required permissions to read and write to specific location in Synapse. Data is stored in the Synapse project at https://www.synapse.org/agora.

To process all Agora data from Synapse files to JSON files suitable for import into the Agora MongoDB, use ([exec/process.R](exec/process.R)), which requires the [config.json](config.json) file. Providing the `--store` parameter will store the results to Synapse.
  
  ```
  exec/./process.R --config config.json --store
  ```

If changes to the gene database used are needed, run [inst/getMyGeneInfo.R](inst/getMyGeneInfo.R) to get gene information from [mygene.info](http://mygene.info). This puts the resulting data into Synapse. This step is required to run separately due to conflicts between the `synapser` and `mygene` R packages.

The output of this tool is used by the [Agora data manager](https://github.com/Sage-Bionetworks/agora-data-manager).

## Development

Create a feature branch from `master` for doing development. You can test the data processing locally using the `exec/process.R` command mentioned above. Note that there are separate configurations for staging and production. You will manually run the processing command when your code is ready, as CICD is not set up in this repository (yet).

### Releases

Use semantic versioning for releases.
