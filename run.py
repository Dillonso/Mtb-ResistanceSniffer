import sys, os
parent_folder = os.path.basename(os.getcwd())
sys.path.append(os.path.join("lib"))
import main

try:
    import psyco
    psyco.profile()
except:
    pass


###############################################################################
if __name__ == "__main__":
    options = {
               "-i":"input",    # input folder
               "-o":"output",   # output folder
               "-f":"",         # parent folder
               "-t":"detailed", # text output detailed | short | no
               "-c":0.75,       # Certainty cutoff in taxonomic splits
               "-g":"",         # generic reference file name
               "-b":"",         # base folder
               "-r":"lib",         # root folder
            }
    options["-g"] = filter(lambda f: f[-4:].upper()==".FAA", os.listdir("sources"))[0][:-4]
    arguments = sys.argv[1:]
    if arguments:
        for i in range(0,len(arguments)-1,2):
            key = arguments[i]
            if key not in options:
                raise IOError("Unknown argument " + key + "!")
            if i <= len(arguments)-2:
                options[key] = arguments[i+1]
    main.execute(options)
