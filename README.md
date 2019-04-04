# Agora Data Processing

Script for processing data from Synapse to Agora.

There are two processes:

1. Get gene information from [mygene.info](http://mygene.info) ([inst/getMyGeneInfo.R](inst/getMyGeneInfo.R))
1. process all Agora data from Synapse files to JSON files suitable for import into the Agora MongoDB ([inst/process.R](inst/process.R) )

## Development

The `master` branch has the most recent changes.

There is a long-running `staging` branch that pushes data to a separate data folder for testing purposes.
