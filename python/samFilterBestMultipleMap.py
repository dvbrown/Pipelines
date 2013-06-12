#sort you SAM file by read ID, so that multiple mappings are in adjacent lines and the write a script to filter the best one
#Written by Simon Anders
import sys, re
import HTSeq

insam = HTSeq.SAM_Reader( sys.stdin )

# Go through all reads, with their alignments bundled up:
for bundle in HTSeq.bundle_multiple_alignments( insam ):
   bestAlmt = None
   # Go through all alignments of a given read, looking
   # for the one with the best alignment score
   for almt in bundle:
      if bestAlmt is None:
         bestAlmt = almt
      elif almt.aQual > bestAlmt.aQual:
         bestAlmt = almt
      elif almt.aQual == bestAlmt:
         # If there are more than one best alignment, 
         # better skip the read
         bestAlmt = None
   if bestAlmt is not None:
      # Change the NH field to 1 and print the line
      print re.sub( "NH:i:\d+", "NH:i:1", bestAlmt.original_sam_line )
      
#call this script with the command sort samfile.sam | python chooseBest.py > filtered.sam	