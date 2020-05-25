#input arguments: sequence

library(NuPoP)
args<-commandArgs(TRUE)

#save the input sequence as a temporary fasta file
sequence = args[1]
fasta = paste(">seq",sequence,sep = '\n')
write(fasta,file = "nupop_temp_files/temp_nucleo_input.fa")

#first compute nucleosome occupancy
capture.output(predNuPoP("nupop_temp_files/temp_nucleo_input.fa",species=7,model=4),file='NUL')

#read the results
results = read.table("temp_nucleo_input.fa_Prediction4.txt", header = TRUE)

#delete the temporary files
invisible(file.remove("NUL","nupop_temp_files/temp_nucleo_input.fa", "temp_nucleo_input.fa_Prediction4.txt"))

#return the results
cat(results[,"Occup"])
