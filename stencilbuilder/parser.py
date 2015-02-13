#!/usr/bin/python
import sys
import shutil
import StringIO

from StencilBuilder import *

# Retrieve the file from the command line arguments.
argc = len(sys.argv)
if (not 2 <= argc < 4):
  raise RuntimeError("Illegal number of arguments")
else:
  inputfile = str(sys.argv[1])
  if (argc == 3):
    outputfile = str(sys.argv[2])

# Save a backup of the original file.
shutil.copyfile(inputfile, "{0}{1}".format(inputfile, ".orig"))

# Read the file into the memory.
f = file(inputfile, "r")
lines = f.readlines()
f.close()

# Remove the line breaks.
for n in range(len(lines)):
  lines[n] = lines[n].rstrip()

record = False
blocks = []

lineindex = -1
block = []
names = []

# Loop over lines and store blocks of StencilBuilder code.
for n in lines:
  lineindex += 1
  line = n.find("//$")

  # Process the keyword
  if(line != -1):
    indent = line
    kw = n[line + 3::].strip().split(' ')
    if(kw[0] == "SBStart"):
      startline = lineindex
      record = True
      continue

    elif(kw[0] == "SBEnd"):
      endline = lineindex
      record = False
      blocks.append((block, startline, endline, indent))
      block = []
      names = []
      continue

    else:
      raise RuntimeError("Syntax error in StencilBuilder tags")

  if(record):
    block.append(n.strip()+'\n')

# Execute the blocks in reverse order
blocks.reverse()

for n in blocks:
  # Get the code block.
  code = ''.join(n[0])

  # Execute the code block and capture the output by temporarily
  # redirecting the stdout.
  buffer = StringIO.StringIO()
  sys.stdout = buffer
  exec(code)
  sys.stdout = sys.__stdout__
  output = buffer.getvalue().splitlines()
  del(code)
  del(buffer)

  # Add the indentation.
  for l in range(len(output)):
    output[l] = '{0}{1}'.format(' ' * n[3], output[l])

  # Delete the StencilBuilder lines.
  del(lines[n[1]:n[2]+1])

  # Replace it with the new code.
  lines[n[1]:n[1]] = output

if (argc == 2):
  for n in lines:
    print("{0}".format(n))
elif (argc == 3):
  f = file("{0}".format(outputfile), "w")
  for n in lines:
    f.write("{0}\n".format(n))
  f.close()
