# runs all the matlab scripts

from subprocess import Popen



def main():
  numbers = [3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,1000]
  argument = ["matlab", "-nodesktop", "-nodisplay", "-nojvm", "-nosplash"]
  for i in numbers:
    fileName = "command" + str(i) + ".txt"
    file = open(fileName)
    Popen(argument, stdin=file)
  return 0

main()
