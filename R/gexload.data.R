gexload.data <-
function (source, genome)
{
   ###############################################################################
   #                                                                             #
   #                                LOAD.DATA                                    #
   #                                                                             #
   #      ACTION = Load the data.Rdata file                                      #
   #      OUTPUT = Matrix "data" of all ENSEMBL annotations                      #
   #                                                                             #
   ###############################################################################

    write("Loading localization data", file = "")
    data.Rdata.dir = file.path(paste(source, genome, ".Rdata", sep = ""))
    if (!file.exists(data.Rdata.dir))
    {
        write(sprintf("ERROR: file ", data.Rdata.dir, " does not exist in the default directory", source), file = "")
        write(paste("Load manually the ", genome, ".Rdata file", sep = ""), file = "")
        load(file = file.choose())
    } else
    {
        load(file = data.Rdata.dir)
        write(paste(data.Rdata.dir, " file Loaded", sep = ""), file = "")
    }
    return(data)
}

