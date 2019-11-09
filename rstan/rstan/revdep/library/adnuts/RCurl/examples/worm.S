counter =
function(max = 10)
{
  counts = numeric()
  numRecords = 0
  
  read = function(txt) {
     
     id = gsub("^>(.*) .*", "\\1", txt)

     vals = table(id)
     
     for(id in names(vals))
       counts[id] <<- ifelse(id %in% names(counts), counts[id], 0) + vals[id]

     numRecords <<- numRecords + 1

     if(numRecords >= max) {
       msg = list(message = "Reached maximum count")
       class(msg) = c("CountException", "error", "condition")
       stop(msg)
     }

#     sum(nchar(txt))
  }

  list(read = read, counts = function() counts)
}  

theCount = counter()
f = chunk2LineReader(theCount$read)

x = postForm("http://www.wormbase.org/db/searches/advanced/dumper",
         species = "elegans", 
         list = "AC3.8  M7  X  IV  III:1000..4000",
         feature = "BLASTX Hits",
         DNA = "feature coordinates only",
         flank5 = "0",
         flank3 = "0",
         relative = "Chromosome",
         orientation = "Relative to feature",
         dump  = "Plain TEXT",
         DUMP = "DUMP",
         .cgifields = c("feature", "orientation", "dump", "relative", "DNA"),
         .opts = list(writefunction = f$read))

nchar(x)





#         DNA = 'feature and flanking sequences',
#         DUMP = 'DUMP',
         #species = 'elegans',

x = postForm('http://www.wormbase.org/db/searches/advanced/dumper',
         species="briggsae",
         list="AC3.8 M7 X IV III:1000..4000",
         flank3="0",
         flank5="0",
         feature="Gene Models",
         dump = "Plain TEXT",
         orientation = "Relative to feature",
         relative = "Chromsome",
         DNA ="flanking sequences only",
         .cgifields = paste(c("feature", "orientation", "DNA", "dump", "relative"), collapse=", "))

print(x)
