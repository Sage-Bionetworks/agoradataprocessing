# Agora Data Processing

Script for processing data from Synapse to Agora. This required permissions to read and write to specific location in Synapse.

To process all Agora data from Synapse files to JSON files suitable for import into the Agora MongoDB, use ([exec/process.R](exec/process.R)), which requires the [config.json](config.json) file. Providing the `--store` parameter will store the results to Synapse.
  
  ```
  exec/./process.R --config config.json --store
  ```

If changes to the gene database used are needed, run [inst/getMyGeneInfo.R](inst/getMyGeneInfo.R) to get gene information from [mygene.info](http://mygene.info). This puts the resulting data into Synapse. This step is required to run separately due to conflicts between the `synapser` and `mygene` R packages.

## Development

The `master` branch has the most recent changes used in production.

The `develop` branch has the most recent changes that are in testing.

There is a long-running `staging` branch that pushes data to a separate data folder for testing purposes.
