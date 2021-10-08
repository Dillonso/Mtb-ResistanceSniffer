import os, string, seq_io, math
import collection, blast, progressbar

global genbank_extensions
genbank_extensions = ['.GB', '.GBK', '.GBF', '.GBFF']
global consensus_extensions
consensus_extensions = ['.CONS', '.CONSENS', '.CONSENSUS']
global contig_extensions
contig_extensions = ['.FASTA', '.FAS', '.FNA', '.FST', ".FA"]
global protein_extensions
protein_extensions = ['.FAA']
global gene_extensions
gene_extensions = ['.FNN']
global vcf_extensions
vcf_extensions = ['.VCF']
global oIO
oIO = seq_io.IO()

# PARA
global FLG_HEURISTIC
FLG_HEURISTIC = False
global PRINT_GRAPH
PRINT_GRAPH = True
global GRAPH_TYPE
GRAPH_TYPE = "HISTOGRAM"

# Thresholds
global e_threshold
e_threshold=0
global score_threshold
score_threshold=0
global mismatch_threshold
mismatch_threshold=0
global coverage_threshold
coverage_threshold=0
global foldchange_threshold
foldchange_threshold=0

global format_number
def format_number(n,dec=2):
    return str(int(n*10**dec)/float(10**dec))

def execute(options):
    def get_para():
        para = {"heuristic":"Yes",
                "graph":"histogram",
                "size":"1"}
        inf_path = os.path.join(options['-b'],options['-f'],"sources","para.inf")
        if os.path.exists(inf_path):
            inf = map(lambda item: item.split("="), oIO.open_text_file(inf_path,True,True))
            for item in inf:
                if item[0] in para:
                    para[item[0]] = item[1]
        return para
        
    para = get_para()
    process(para,options)    

def process(para,options):
    def check_type(fname):
        try:
            extension = fname[fname.rfind("."):].upper()
        except:
            return
        if extension in genbank_extensions:
            return "genbank"
        elif extension in consensus_extensions:
            return "consensus"
        elif extension in contig_extensions:
            return "contig"
        elif extension in protein_extensions:
            return "protein"
        elif extension in gene_extensions:
            return "gene"
        elif extension in vcf_extensions:
            return "VCF"
        else:
            return
        
    def parse_table(path=""):
        path = os.path.join(path,"sources","table.txt")
        f = open(path)
        data = f.read().split("\n")
        f.close()
        oTable = None
        index = []
        first = ""
        abbreviation_mode = statistics_mode = False
        for line in data:
            line = line.strip()
            if not line or line[0] in ("'","\"","#"):
                continue
            if line.strip() == "//":
                abbreviation_mode = statistics_mode = False
                continue
            if line.strip() == "// Abbreviations":
                abbreviation_mode = True
                continue
            if line.strip() == "// Statistics":
                statistics_mode = True
                continue
            
            if abbreviation_mode:
                oTable.add_abbreviation(map(lambda s: s.strip(),line.split("\t")))
                continue
            if statistics_mode:
                oTable.add_statistics(map(lambda s: s.strip(),line.split("\t")))
                continue
            if line.find("Split") > -1:
                line = line.split("\t")
                index = line[0][6:-1].split(".")[1:]
                first,second = line[1].split("|")
                description = ""
                if len(line) > 2:
                    description = line[2]
                if oTable == None:
                    oTable = identification_key("Root",1,first,second)
                else:
                    oTable.add_child_nodes(index,first,second,description)
                continue
            try:
                id,position,power,alleles,values,locus,gene,codon,annotation = line.split("\t")[:9]
                oTable.add_entry(index,key_entry(int(id),int(position),float(power),alleles,values,locus,codon))
            except:
                oTable.add_entry(index,key_entry("00")) #   add empty entry
        return oTable
    
    def parse_reference(path,generic_name,flg_dna=False):
        seq,p = oIO.openFasta(os.path.join(path,"sources",generic_name+".fasta"))
        REFGENOME_LENGTH = len(seq.values()[0])
        if flg_dna:
            path = os.path.join(path,"sources",generic_name+".fnn")
        else:
            path = os.path.join(path,"sources",generic_name+".faa")
        refseq,path = oIO.openFasta(path)
        return refseq
    
    global FLG_HEURISTIC
    global PRINT_GRAPH
    global GRAPH_TYPE
    global REFGENOME_LENGTH
    global CERTAINTY_CUTOFF
    CERTAINTY_CUTOFF = options["-c"]
    if para["heuristic"] == "Yes":
        FLG_HEURISTIC = True
    else:
        FLG_HEURISTIC = False
    if para["graph"] and para["graph"].upper() != "NO":
        PRINT_GRAPH = True
    else:
        PRINT_GRAPH = False
    GRAPH_TYPE = para["graph"].upper()
    
    #### MAIN
    # from upper input folder
    infolder = os.path.join(options['-f'],options['-i'])
    # from current input folder
    #infolder = options['-i']
    
    outfolder = os.path.join(options['-b'],options['-f'],options['-o'])
    fnames = filter(lambda fname: check_type(fname),os.listdir(os.path.join(options['-b'],infolder)))
    
    if not fnames:
        # Try to read from the general input folder
        infolder = os.path.join(options['-b'],options['-i'])
        fnames = filter(lambda fname: check_type(fname),os.listdir(infolder))
    if not fnames:
        return
    oTable = parse_table(os.path.join(options['-b'],options['-f']))
    
    winners = []
    members = []
    genome_counter = 0
    file_type = ""
    if para['graph'] and para['graph'].upper() != "NO":
        svg = get_svg_template(options,para['graph'])
    for fname in fnames:
        print
        print fname
        oMessenger = None
        new_file_type = check_type(fname)
        if new_file_type != file_type:
            file_type = new_file_type
            flg_dna = False
            if file_type in ("contig","VCF"):
                flg_dna = True 
            oReference = parse_reference(os.path.join(options['-b'],options['-f']),options['-g'],flg_dna)
        obj = None
        if file_type == "genbank":
            obj = Genbank(options['-r'],options['-b'],os.path.join(infolder,fname),oTable,oReference)
        elif file_type == "consensus":
            obj = Consensus(options['-r'],options['-b'],os.path.join(infolder,fname),oTable,oReference)
        elif file_type == "contig":
            obj = Contigs(options['-r'],options['-b'],os.path.join(infolder,fname),oTable,oReference)
        elif file_type == "protein":
            obj = Protein(options['-r'],options['-b'],os.path.join(infolder,fname),oTable,oReference)
        elif file_type == "gene":
            obj = Gene(options['-r'],options['-b'],os.path.join(infolder,fname),oTable,oReference)
        elif file_type == "VCF":
            obj = VCF(options['-r'],options['-b'],os.path.join(infolder,fname),oTable,oReference)
        oMessenger = obj()
        genome_counter += 1
        if oMessenger != None:
            if options["-t"].upper() != "NO":
                flg_detailed = True
                if options["-t"].upper() == "SHORT":
                    flg_detailed = False
                report = oMessenger.tostring(flg_detailed)
                oIO.save(report,os.path.join(options["-b"],options["-f"],"output",fname[:fname.rfind(".")]+".txt"))
            if para['graph'].upper() == "CLADOGRAM": 
                if genome_counter%int(para['size'])==0:
                    if int(para['size']) == 1:
                        svg = cladogram_svg(svg,[fname],[oMessenger.get_winner()],report,True)
                    else:
                        winners.append(oMessenger.get_winner())
                        members.append(fname)
                        svg = add_graphs_to_svg(svg,members,winners,report)
                        winners = []
                        members = []
                    if int(para['size']):
                        save_svg(svg,fname,int(para['graph']),genome_counter,os.path.join(options["-b"],options["-f"]))
                        svg = get_svg_template(options,int(para['graph']))
                else:
                    winners.append(oMessenger.get_winner())
                    members.append(fname)
            elif para['graph'].upper() == "HISTOGRAM":
                svg = histogram_svg(get_svg_template(options,para['graph']),"%s: %s" % (fname,oMessenger.get_top_scored_clade()),oTable,
                                    oMessenger.get_properties(),oMessenger.default_property,
                                    [oMessenger.default_property] + oMessenger.properties)
                save_svg(svg,fname,int(para['size']),genome_counter,os.path.join(options["-b"],options["-f"]))
               
    if para['graph'] and para['graph'].upper() != "NO" and len(members) and FLG_HEURISTIC:
        svg = cladogram_svg(svg,members,winners,report)
        save_svg(svg,fname,int(para['graph']),genome_counter,os.path.join(options["-b"],options["-f"]))
    
def get_svg_template(options,graph):
    if graph.upper() == "CLADOGRAM":
        return oIO.open_text_file(os.path.join(options['-b'],options['-f'],"sources","cladogram.svg"),True)
    elif graph.upper() == "HISTOGRAM":
        return oIO.open_text_file(os.path.join(options['-b'],options['-f'],"sources","histogram.svg"))
    return ""

def save_svg(svg,fname,number_of_graphs,number_of_genomes,basefolder=""):
    svg_fname = fname[:fname.rfind(".")]+".svg"
    if number_of_graphs > 1:
        start = number_of_genomes-number_of_graphs+1
        if start <= 0:
            start = 1
        svg_fname = "plot_%d-%d.svg" % (start,number_of_genomes)
    oIO.save("\n".join(svg),os.path.join(basefolder,"output",svg_fname))
        
def cladogram_svg(svg,members,winners,report,flg_singular=False):
    def get_coord(clade,count,clades):
        x_shift = 0
        y_shift = 10
        y = clades[clade][1]
        if count > 2:
            y_shift = 8+float(clades[clade][2])
            x_shift = 10+float(clades[clade][2])
            count -= 3
        x = x_shift + clades[clade][0]+2.0*float(clades[clade][2])*count
        y += y_shift
        return x,y
    def get_clades(svg):
        height = int(str(svg[0][svg[0].rfind(" "):-2]).strip('"'))
        svg = "\n".join(svg)
        block = map(lambda item: item.split(","), svg[svg.find("<!--"):svg.find("-->")-1].split("\n")[1:])
        return height,dict(zip(map(lambda item: item[0], block),map(lambda item: map(lambda v: int(v), item[1:])+[9], block)))
    def parse_report(report):
        block = map(lambda item: item.split("\t"), report[report.find("Alleles:"):report.find(">Barcode sequence")].split("\n"))[2:-1]
        count = left = right = matches = conflicts = mismatches = 0
        for line in block:
            if line[0]:
                matches += max([left,right])
                conflicts += min([left,right])
                left = right = 0
                continue
            if len(line) < 2 or not line[1]:
                continue
            count += 1
            l,r = map(lambda v: float(v), line[4][:-1].split("/"))
            left += l
            right += r
            if not l and not r:
                mismatches += 1
        matches += max([left,right])
        conflicts += min([left,right])
        return count,matches,conflicts,mismatches
    
    height,clades = get_clades(svg)
    for i in range(len(winners)):
        clade = winners[i]
        
        if clade[:12]=="Not detected":
            continue
        
        x,y = get_coord(clade,winners.count(clade)-winners[i:].count(clade),clades)
        if not x:
            continue
        total,matches,conflicts,mismatches = parse_report(report)
        if not total:
            continue
        if not flg_singular:
            svg.insert(-1,("<circle fill=\"#F7F619\" stroke=\"#000000\" title=\"%s\" cx=\"%f\" cy=\"%d\" r=\"8\"/>\n" % (members[i],x,y))+
                    "<text x=\"%f\" y=\"%d\" font-size=\"%d\" text-anchor=\"middle\">%d</text>" % (x,y+2,10,i+1)+
                    "<text x=\"%f\" y=\"%d\" font-size=\"%d\" text-anchor=\"start\">%d. %s     -     clade %s; matches - %d (%d %%), conflicts - %d (%d %%), mismatches - %d (%d %%)</text>" % 
                    (15,height+20*(i+4),16,i+1,members[i],clade,matches,100-int(100.0*(conflicts+mismatches)/total),
                    conflicts,int(100.0*conflicts/total),mismatches,int(100.0*mismatches/total)))
        else:
            svg.insert(-1,("<circle fill=\"#F7F619\" stroke=\"#000000\" title=\"%s\" cx=\"%f\" cy=\"%d\" r=\"8\"/>\n" % (members[i],x,y)))
    svg[0] = svg[0][:svg[0].rfind(" ")] + (" %s" % (str(height+20*(i+4)))) + svg[0][-2:]
    return svg

def histogram_svg(svg,title,oTable,properties=[],default_property="",property_ls=[],font_size=14):   # properties = [[key,value],...]
    def add_whisker(X,Y,halfheight):
        ls = []
        y1 = Y-halfheight
        if y1 < 51:
            y1 = 51
        y2 = Y+halfheight
        if y2 > 498:
            y2 = 498
        ls.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"%s\" stroke-width=\"%d\" />" % 
                (X,y1,X,y2,"black",2))
        return ls
    
    svg = svg.split("\n")[:-3]
    if not properties:
        properties = {'Sensitive': [1.0,0.0]}
        #svg.append("</svg>")
        #return svg
    if "Sens" in properties:
        properties["Sensitive"] = [properties["Sens"][0],properties["Sens"][1]]
        del properties["Sens"]
    if not property_ls:
        property_ls = properties.keys()
    #properties = dict(zip(map(lambda item: item[0], properties),map(lambda item: [item[1],item[2]], properties)))
    segment = 600/len(property_ls)
    width = segment/2
    svg.append(("<text x=\"50\" y=\"50\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"white\" font-size=\"%d\" " +
        "style=\"text-anchor:start\">%s</text>") %
        (font_size+2,title))
    for i in range(len(property_ls)):
        tag = property_ls[i]
        height = 0
        color = "red"
        description = oTable.get_abbreviation(tag)
        if tag in properties:
            if properties[tag][0] == "ND":
                height = 20
                color = "red"
                description += ": data insufficient"
            elif properties[tag][0] == "DEFAULT":
                height = 20
                color = "grey"
                description += ": by default"
            else:
                if properties[tag][0] < 0.7:
                    color = "orange"
                if properties[tag][0] < 0.3:
                    color = "white"
                height = properties[tag][0]*400
        x = 100 + segment/2 + i*segment
        if tag == default_property:
            try:
                description += ": %s%%" % format_number(oTable.get_statistics(properties[tag][0]),1)
            except:
                pass
            color = "lightgreen"
            svg.append("<a xlink:title=\"%s\"><rect x=\"%d\" y=\"%d\" fill=\"%s\" strock=\"black\" width=\"%d\" height=\"%d\" /></a>" % 
                (description,100 + segment/4 + i*segment,499-height,color,width,height))
            svg.append(("<text x=\"%d\" y=\"532\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"white\" font-size=\"%d\" " +
                "style=\"text-anchor:middle\" transform=\"rotate(90 %d,532)\">%s</text>") %
                (x-3,font_size,x-3,tag))
        else:
            try:
                description += ": resistance %s%%" % format_number(oTable.get_statistics(properties[tag][0]),1)
            except:
                pass
            svg.append("<a xlink:title=\"%s\"><rect x=\"%d\" y=\"%d\" fill=\"%s\" strock=\"black\" width=\"%d\" height=\"%d\" /></a>" % 
                (description,100 + segment/4 + i*segment,499-height,color,width,height))
            svg.append(("<text x=\"%d\" y=\"520\" font-family=\"Times New Roman\" font-weight=\"bold\" fill=\"white\" font-size=\"%d\" " +
                "style=\"text-anchor:middle\">%s</text>") %
                (x,font_size,tag))
        if properties[tag][1]:
            svg += add_whisker(100 + segment/2 + i*segment,499-height,properties[tag][1]*200)
        
    svg.append("</svg>")
    return svg

###############################################################################
class identification_key(collection.Collection):
    def __init__(self,name,index,first_child=None,second_child=None):
        collection.Collection.__init__(self,name,index)
        self.child_nodes = []
        self.descriptions = []
        self.power_values = []
        self.abbreviations = {}
        self.statistics = []
        self._set_child_nodes(first_child,second_child)
        
    def _do(self,request="",arguments=[]):
        if request == "is resistance unit":
            return self.is_resistance_unit()
        elif request == "get title":
            return self.title
        elif request == "get child nodes":
            return self.child_nodes
        elif request == "get child node titles":
            return self.get_child_node_titles()
        else:
            pass
            
    def _set_child_nodes(self,first_child,second_child,description=""):
        if first_child == None:
            return
        self.child_nodes.append(key_unit(first_child,len(self.child_nodes)+1))
        self.child_nodes.append(key_unit(second_child,len(self.child_nodes)+1))
        if description:
            data = description.split("=")
            self.descriptions.append(data[0])
            try:
                self.power_values.append(float(data[1]))
            except:
                self.power_values.append(1.0)
        else:
            self.descriptions.append("")
            self.power_values.append(1.0)
    
    def add_child_nodes(self,index,first_child=None,second_child=None,description=""):
        if index:
            i = int(index[0])
            index = index[1:]
            self.child_nodes[i-1].add_child_nodes(index,first_child,second_child,description)
            return
        self._set_child_nodes(first_child,second_child,description)
        
    def add_entry(self,index,oEntry):
        if index:
            i = int(index[0])
            index = index[1:]
            self.child_nodes[i-1].add_entry(index,oEntry)
            return
        oEntry.pair = len(self.child_nodes)/2
        oEntry.trigger = self._do
        self.append(oEntry)
        
    def add_abbreviation(self,item):
        try:
            abbr,term = item
        except:
            return
        self.abbreviations[abbr] = term
        
    def add_statistics(self,item):
        try:
            rank,val = map(lambda v: float(v),item)
        except:
            return
        self.statistics.append([rank,val])
        if len(self.statistics) > 1:
            self.statistics.sort()
            
    def get_abbreviation(self,key):
        if key in self.abbreviations:
            return self.abbreviations[key]
        return key
    
    def get_statistics(self,val,flg_reverse=False):
        try:
            val = float(val)
        except:
            return
        for i in range(1,len(self.statistics),1):
            if val < self.statistics[i][0]:
                if flg_reverse:
                    return 100.0 - self.statistics[i-1][1]
                return self.statistics[i-1][1]
        if flg_reverse:
            return 100.0 - self.statistics[-1][1]
        return self.statistics[-1][1]
                
    def get_child_node_titles(self):
        if not self.child_nodes:
            return ""
        return " | ".join(map(lambda child_node: child_node.title, self.child_nodes))
    
    def get_power(self):
        return self.power
        
    def trace_path(self,index=1):
        if not len(self):
            return self.title
        return "%s | %s" % (self.title,self.child_nodes[index-1].trace_path(index))
    
    def is_resistance_unit(self):
        return len(filter(lambda title: title.find("~") > 0, self.get_child_node_titles().split(" | ")))
        
    def identify(self,mode,basefolder,dbname,source_path,query_fasta,oMessenger=None):
        def identify_by_protalignment(oMessengerUnit,basefolder,dbname,source_path,query_fasta):
            fasta,path = oIO.openFasta(source_path,True)
            bar = progressbar.indicator(len(self),"BLAST")
            counter = 1
            self.container.sort()
            try:
                for k in range(len(self)):
                    entry = self.container[k]
                    
                    if entry.is_empty():
                        oMessengerUnit.clade = oMessengerUnit.title
                        print "\tIdentified node: " + oMessengerUnit.clade
                        print
                        return oMessengerUnit,True
                    
                    try:
                        query = filter(lambda record: record.split("|")[1].strip() == entry.locus, query_fasta.keys())[0]
                    except:
                        continue
                    oBlast = blast.BLAST("blast","protein",basefolder,dbname,query_fasta[query],os.path.join(basefolder,"tmp"))
                    oBlast.execute()
                    try:
                        e,score,title,hsp = oBlast.get_top_alignment()
                        oHSP = hsp[0]
                    except:
                        continue
                    if (e == None or 
                        (oMessengerUnit.e_threshold and e > oMessengerUnit.e_threshold) or 
                        (oMessengerUnit.score_threshold and score < oMessengerUnit.score_threshold)):
                        continue
                    sbjct = fasta[int(title.split("|")[0])-1]
                    try:
                        oMessengerUnit.append(MessengerEntry(title,entry.copy(),sbjct[entry.codon+oHSP.sbjct_start-oHSP.query_start-1]))
                    except:
                        oMessengerUnit.append(MessengerEntry(title,entry.copy(),"*"))
                    bar(counter)
                    counter += 1
            except:
                bar.stop()
                print
                print "Error in the input file"
                return None,True
            bar.stop()
            return oMessengerUnit,False
        
        def identify_by_dnaalignment(oMessengerUnit,basefolder,dbname,source_path,query_fasta):
            oMessengerUnit.e_threshold = 0.0001
            fasta,path = oIO.openFasta(source_path,True)
            bar = progressbar.indicator(len(self),"BLAST")
            counter = 1
            self.container.sort()
            
            try:
                for k in range(len(self)):
                    entry = self.container[k]
                    if entry.is_empty():
                        oMessengerUnit.clade = oMessengerUnit.title
                        print "\tIdentified node: " + oMessengerUnit.clade
                        print
                        return oMessengerUnit,True
                    
                    try:
                        query = filter(lambda record: record.split("|")[1].strip() == entry.locus, query_fasta.keys())[0]
                    except:
                        continue
                    oBlast = blast.BLAST("blast","dna",basefolder,dbname,query_fasta[query],os.path.join(basefolder,"tmp"))
                    oBlast.execute()
                    try:
                        e,score,title,hsp = oBlast.get_top_alignment()
                        oHSP = hsp[0]
                    except:
                        continue
                    if (e == None or 
                        (oMessengerUnit.e_threshold and e > oMessengerUnit.e_threshold) or 
                        (oMessengerUnit.score_threshold and score < oMessengerUnit.score_threshold)):
                        continue
                    for h in range(len(hsp)):
                        location = (entry.codon-1)*3 + 1
                        if location < hsp[h].query_start or location > hsp[h].query_end:
                            continue
                        oMessengerUnit.add(MessengerEntry(title,entry.copy(),oIO.translate([hsp[h].sbjct[location-hsp[h].query_start:location-hsp[h].query_start+3]])[0]))
                        break
                    bar(counter)
                    counter += 1
            except:
                bar.stop()
                print
                print "Error in the input file"
                return None,True
            bar.stop()
            return oMessengerUnit,False
        
        def identify_by_comparison(oMessengerUnit,basefolder,query_fasta,reference):
            def get_residue(sequence,query,entry,oParser):
                coords = map(lambda v: int(v), query.split("|")[-1].strip()[1:-1].split(".."))
                strand = 'dir'
                if entry.locus[-1] == "c":
                    strand = 'rev'
                    codon_p = (coords[1]-entry.position)%3
                    seq = oParser.translate(sequence[entry.position-3+codon_p:entry.position+codon_p],strand)
                else:
                    codon_p = (entry.position-coords[0])%3
                    seq = oParser.translate(sequence[entry.position-codon_p:entry.position+3-codon_p],strand)
                return seq

            oParser = seq_io.Parser()
            self.container.sort()
            for k in range(len(self)):
                entry = self.container[k]
                if entry.is_empty():
                    oMessengerUnit.clade = oMessengerUnit.title
                    print "\tIdentified node: " + oMessengerUnit.clade
                    print
                    return oMessengerUnit,True
                    
                try:
                    query = filter(lambda record: record.split("|")[1].strip() == entry.locus, reference.keys())[0]
                except:
                    continue
                seq = get_residue(query_fasta,query,entry,oParser)
                oMessengerUnit.append(MessengerEntry(entry.title,entry.copy(),seq))
            return oMessengerUnit,False
            
        def identify_by_polymorphic_sites(oMessengerUnit,basefolder,dataset,reference):
            def get_residue(ref_gene,allele,entry,oParser):
                seq = ref_gene[3]
                codon = seq[3*(entry.codon-1):3*(entry.codon-1)+3]
                p = (entry.position-ref_gene[1]-1)%3
                if ref_gene[0].strip()[-1] == "c":
                    p = abs(p-2)
                    allele = oParser.reverse_complement(allele)
                return oParser.translate(codon[:p]+allele+codon[p+1:])
            
            oParser = seq_io.Parser()
            self.container.sort()
            for k in range(len(self)):
                entry = self.container[k]
                if entry.is_empty():
                    oMessengerUnit.clade = oMessengerUnit.title
                    print "\tIdentified node: " + oMessengerUnit.clade
                    print
                    return oMessengerUnit,True
                    
                ref_gene = filter(lambda record: record[0] == entry.locus, reference)[0]
                seq = ref_gene[3]
                if ref_gene[0].strip()[-1] == "c":
                    seq = oParser.reverse_complement(seq)
                try:
                    ref,allele = filter(lambda record: record[0] == entry.position and record[3] > 20, dataset)[0][1:3]
                except:
                    allele = seq[entry.position-ref_gene[1]-1]
                    ref = allele
                if ref != seq[entry.position-ref_gene[1]-1]:
                    continue
                residue = get_residue(ref_gene,allele,entry,oParser)
                oMessengerUnit.append(MessengerEntry(entry.title,entry.copy(),residue))
            return oMessengerUnit,False

        #    MAIN
        if not len(self):
            return oMessenger
        if oMessenger == None:
            oMessenger = Messenger(mode)
        oMessengerUnit = MessengerUnit(self.title,self.get_child_node_titles(),len(self))
        oMessengerUnit.descriptions = []
        oMessengerUnit.descriptions.extend(self.descriptions)
        oMessengerUnit.power_values = []
        oMessengerUnit.power_values.extend(self.power_values)
        print  "Identify between %s" % oMessengerUnit.format_child_node_titles()
        
        if mode == "CDS":
            oMessengerUnit,flg_finish = identify_by_protalignment(oMessengerUnit,basefolder,dbname,source_path,query_fasta)
        elif mode == "contigs":
            oMessengerUnit,flg_finish = identify_by_dnaalignment(oMessengerUnit,basefolder,dbname,source_path,query_fasta)
        elif mode == "consensus":
            oMessengerUnit,flg_finish = identify_by_comparison(oMessengerUnit,basefolder,query_fasta,dbname)
        elif mode == "vcf":
            oMessengerUnit,flg_finish = identify_by_polymorphic_sites(oMessengerUnit,basefolder,dbname,query_fasta)
        if oMessengerUnit == None:
            return
        if flg_finish:
            oMessenger.append(oMessengerUnit)
            return oMessenger

        if FLG_HEURISTIC:
            oMessenger.append(oMessengerUnit)
            winners = oMessengerUnit.identify()
            if not winners:
                if oMessengerUnit.clade:
                    print "\t" + oMessengerUnit.clade
                return oMessenger
            if len(winners) == 1:
                return self.child_nodes[winners[0]].identify(mode,basefolder,dbname,source_path,query_fasta,oMessenger)
            map(lambda i: oMessenger.child_nodes.append(self.child_nodes[i].identify(mode,basefolder,dbname,source_path,query_fasta,Messenger(mode))),
                winners)
            return oMessenger
        else:
            oMessenger.append(oMessengerUnit)
            map(lambda i: oMessenger.child_nodes.append(self.child_nodes[i].identify(mode,basefolder,dbname,source_path,query_fasta,Messenger(mode))),
                range(len(self.child_nodes)))
            return oMessenger

class key_unit(identification_key):
    def __init__(self,name,index):
        identification_key.__init__(self,name,index)
        
    def copy(self):
        return self
        
class key_entry:
    def __init__(self,id,position=None,power=0,alleles="",values="",locus="",codon=0,pair=1,trigger=None):
        self.id = id
        self.position = None
        try:
            self.position = int(position)
        except:
            pass
        self.power = float(power)
        self.alleles = alleles
        self.values = values
        self.title = self.locus = locus
        self.codon = int(codon)
        self.pair = pair
        self.trigger = trigger
        
    def __eq__(self,other):
        if type(other) != type(self):
            return
        if self.id == other.id and self.alleles == other.alleles:
            return True
        return
    
    def __cmp__(self,other):
        if type(other) != type(self):
            return cmp(self,other)
        if self.id < other.id:
            return -1
        if self.id > other.id:
            return 1
        if self.id == other.id and self.position == other.position and self.alleles == other.alleles:
            return 0
        else:
            return cmp(self.alleles,other.alleles)
        
    def get_alleles(self):
        return map(lambda item: item.split(","),self.alleles.split("|"))
    
    def get_values(self):
        alleles = self.get_alleles()
        values = map(lambda item: item.split(","),self.values.split("|"))
        return map(lambda i: dict(zip(alleles[i],map(lambda v: float(v), values[i]))),range(len(alleles)))
    
    def get_score(self,i,L):
        values = self.get_values()
        if L not in values[i]:
            return 0
        return values[i][L]
        
    def get_maxScore(self,allele):
        if allele in reduce(lambda a,b: a+b, map(lambda item: item.split(","), self.alleles.split("|"))):
            return max(reduce(lambda a,b: a+b, map(lambda item: map(lambda v: float(v), item.split(",")), self.values.split("|"))))
        return 0
    
    def get_power(self,i):
        i -= 2*(self.pair-1)
        values = self.get_values()
        return max(map(lambda L: values[i][L], values[i].keys()))
    
    def is_empty(self):
        if self.position == None:
            return True
        return False
    
    def copy(self):
        return key_entry(self.id,self.position,self.power,self.alleles,self.values,self.locus,self.codon,self.pair,self.trigger)
        
class Messenger(collection.Collection):
    def __init__(self,input_type="",title=""):
        collection.Collection.__init__(self,title)
        self.input_type = input_type
        self.e_threshold = e_threshold
        self.score_threshold = score_threshold
        self.success_cutoff = 5
        self.child_nodes = []
        self.top_scorred_clade = None
        self.properties = ["EMB","FLQ","INH","RIF","PAS","SM","KAN","CM","PZA","AMK","CS","ETA","OFL"]
        self.properties.sort()
        self.default_property = "Sensitive"
       
    def __add__(self,other):
        pass

    def get_properties(self,default_property="",default_value=1.0):
        if self.top_scorred_clade == None:
            self.top_scorred_clade = self.get_top_scored_clade()
        properties = self.top_scorred_clade.get_predicted_resistance(self.properties)
        if default_property:
            self.default_property = default_property
        return properties
        
    def get_lineage_path(self):
        lineage = map(lambda i: self[i].clade.split(" | ")[0], range(len(self)))
        for i in range(len(lineage)-1,0,-1):
            if lineage[i]==lineage[i-1] or lineage[i].find("~") > -1:
                del lineage[i]
        return " | ".join(lineage)
    
    def get_polymorphic_sites(self):
        output = []
        for oMU in self:
            for oSubMU in oMU.split_entries():
                output.append("Node %s, distinguishing between %s %s" % (oSubMU.title,oSubMU.child_node_titles,oSubMU.descriptions[0]))
                output.append("\t"+"\n\t".join(oSubMU.enroll()))
        if self.child_nodes:
            child_nodes = oMU.get_childnode_titles()
            output.append("SPLIT BETWEEN %s" % " AND ".join(child_nodes))
            for i in range(len(self.child_nodes)):
                output.append("\tSUBTREE %s\n" % child_nodes[i])
                output += self.child_nodes[i].get_polymorphic_sites()
                output.append("#"*20)
        return output

    def tostring(self,flg_detailed="True"):
        output=[]
        if not len(self):
            return "Identification failed!"

        if self.top_scorred_clade == None:
            self.top_scorred_clade = self.get_top_scored_clade()
        output.append("Clade(s): " + self.top_scorred_clade.tostring(self.properties))

        output.append("\nLineages:")
        output.append(self.get_lineage_path())
        if self.child_nodes:
            p = len(output[-1])
            for oME in self.child_nodes:
                output[-1] += " | " + oME.get_lineage_path() + ("\n" + "."*p)
            output[-1] = output[-1].strip(".")
                
        if not flg_detailed:
            return "\n".join(output)

        output.append("\n\nAlleles:\n") 
        output += self.get_polymorphic_sites()
        
        return "\n".join(output)
    
    def get_score(self,counter=0,score=0):
        if not len(self):
            return []
        titles=[]
        for oMU in self:
            if oMU.is_resistance_unit():
                score /= counter
                return [ScoreUnit(oMU.title,score,oMU.get_found_sites(),oMU.get_sensitivity_pannel(),oMU.get_resistance_pannel())]
            if oMU.is_empty():
                score /= counter
                return [ScoreUnit(oMU.title,score,oMU.get_found_sites(),oMU.get_sensitivity_pannel(),oMU.get_resistance_pannel())]
            score += max(oMU.scores)
            titles = oMU.get_childnode_titles()
            counter += 1
        if self.child_nodes:
            return reduce(lambda a,b: a+b, map(lambda i: self.child_nodes[i].get_score(counter,score), 
                    range(len(self.child_nodes))))
        #print oMU.scores,counter,titles,oMU.title,titles.index(oMU.title.strip())
        #scores = map(lambda v: v/counter, oMU.scores)
        score = oMU.scores[titles.index(oMU.title.strip())]
        return [ScoreUnit(oMU.title.strip(),score)]
    
    def get_top_scored_clade(self):
        if self.top_scorred_clade != None:
            return self.top_scorred_clade
        oScores = Scores()
        oScores.set_scores(self.get_score())
        self.top_scorred_clade = oScores.winner
        return self.top_scorred_clade
    
    def enroll(self,lineage):
        output = []
        if not lineage:
            return output
        output.append("Node %s, distinguishing between %s" % (self[0].title,self[0].child_node_titles))
        i = self[0].child_node_titles.split(" | ").index(lineage[0])
        scores = [0,0]
        for entry in self[0]:
            allele = entry.allele
            alleles = entry.entry.alleles.split("|")[i]
            if allele not in alleles.split(","):
                allele = allele.lower()
            output.append("\t%d\t%s [%d]\t%s [%s]\t%s/%s;" % (int(entry.entry.position),
                entry.entry.locus,int(entry.entry.codon),allele,alleles,
                format_number(entry.entry.get_score(0,allele)),
                format_number(entry.entry.get_score(1,allele))))
            scores[0] += entry.entry.get_score(0,allele)
            scores[1] += entry.entry.get_score(1,allele)
        output.append("\t\t\t\t\t\t\t%s/%s" % (format_number(scores[0]/len(self[0])),format_number(scores[1]/len(self[0]))))
        output.append("")
        return output + self.child_nodes[i].enroll(lineage[1:])

    def get_winner(self,clade=""):
        if FLG_HEURISTIC and len(self):
            return self[-1].clade.split("~")[0]
        if not self.child_nodes:
            return clade
        i = self[0].scores.index(max(self[0].scores))
        return self.child_nodes[i].get_winner(self[0].clade)
    
    def get_sequence(self,lineage):
        if not lineage:
            return ""
        seq = self[0].get_sequence()
        i = self[0].child_node_titles.split(" | ").index(lineage[0])
        return seq + self.child_nodes[i].get_sequence(lineage[1:])
    
    def copy(self):
        oMessenger = Messenger(self.input_type,self.title)
        oMessenger.e_threshold = self.e_threshold
        oMessenger.score_threshold = self.score_threshold
        oMessenger.success_cutoff = self.success_cutoff
        oMessenger.child_nodes = map(lambda oMU: oMU.copy(), self.child_nodes)
        return oMessenger

class MessengerUnit(Messenger):
    def __init__(self,title,child_node_titles,size):
        Messenger.__init__(self,"",title)
        self.size = size
        self.child_node_titles = child_node_titles
        self.mismatch_threshold = mismatch_threshold
        self.coverage_threshold = coverage_threshold
        self.foldchange_threshold = foldchange_threshold
        self.scores = []
        self.descriptions = []
        self.power_values = []
        self.clade = ""
        
    def add(self,oMessengerEntry):
        self.append(oMessengerEntry)
        
    def identify(self):
        if not len(self) or not self.size or not self.child_node_titles or (self.coverage_threshold and self.coverage_threshold > float(len(self))/self.size):
            return [-1]
        child_nodes = self.child_node_titles.split(" | ")
        if len(child_nodes) == 1:
            return [0]
        self.scores = map(lambda i: 0, range(len(child_nodes)))
        mismatches = 0
        scored_power = 0
        if not len(self):
            return [-1]
        for oME in self:
            alleles = oME.entry.get_alleles()
            success = False
            for i in range(len(alleles)):
                if oME.allele.upper() in alleles[i]:
                    success = True
                    self.scores[i+2*(oME.entry.pair-1)] += oME.entry.get_score(i,oME.allele.upper())
            if not success:
                mismatches += 1
        self.scores = map(lambda i: self.normalize_score(i),range(len(self.scores)))
        if self.mismatch_threshold and mismatches > self.mismatch_threshold:
            self.clade = "Not detected | Mismatch threshold %d > %d" % (mismatches,self.mismatch_threshold)
            return [-1]
        if min(self.scores)==max(self.scores):
            self.clade = "Not detected | Fold change threshold %f < %f" % (1.00,self.foldchange_threshold)
            return [-1]
        if self.foldchange_threshold and min(self.scores) and float(max(self.scores))/min(self.scores) < self.foldchange_threshold:
            self.clade = "Not detected | Fold change threshold %f < %f" % (float(max(self.scores))/min(self.scores),self.foldchange_threshold)
            return [-1]
        ordered_scores = self.get_ordered_scores()
        i = self.scores.index(ordered_scores[0])
        self.clade = child_nodes[i]
        if len(ordered_scores) < 2 or ordered_scores[0] > CERTAINTY_CUTOFF or self[0].entry.trigger("is resistance unit"):
            return [i]
        j = self.scores.index(ordered_scores[1])
        return [i,j]
    
    def get_score(self):
        if not len(self):
            return 0
        return sum(map(lambda oMessengerEntry: oMessengerEntry.get_maxScore(), self))
    
    def get_ordered_scores(self):
        if not len(self):
            return []
        scores = map(lambda v: v, self.scores)
        scores.sort()
        scores.reverse()
        return scores
    
    def is_discriminative_power_sufficient(self):
        return self.get_score() >= self.success_cutoff
    
    def is_resistance_unit(self):
        return len(filter(lambda title: title.find("~") > 0, self.get_childnode_titles()))
    
    def is_empty(self):
        if len(self)+len(self.scores):
            return False
        return True
    
    def get_found_sites(self):
        if not self.is_resistance_unit():
            return {}
        titles = self.get_childnode_titles()
        if len(self.scores)==len(titles):
            return dict(zip(
                        map(lambda i: titles[i].split("~")[-1], range(1,len(titles),2)),
                        map(lambda i: len(filter(lambda oME: oME.entry.pair==i+1 and 
                            (oME.entry.get_score(0,oME.allele) or oME.entry.get_score(1,oME.allele)),self)),
                            range(len(self.scores)/2))
                    ))
        # Empty terminal unit with antibiotic resistance set by default
        return dict(zip(
                    map(lambda i: titles[i].split("~")[-1], range(1,len(titles),2)),
                    map(lambda i: 1.0, range(1,len(titles),2))
                ))
    
    def get_resistance_pannel(self):
        if not self.is_resistance_unit():
            return {}
        titles = self.get_childnode_titles()
        if len(self.scores)==len(titles):
            return dict(zip(
                        map(lambda i: titles[i].split("~")[-1], range(1,len(titles),2)),
                        map(lambda i: self.scores[i], range(1,len(self.scores),2))
                    ))
        # Empty terminal unit with antibiotic resistance set by default
        return dict(zip(
                    map(lambda i: titles[i].split("~")[-1], range(1,len(titles),2)),
                    map(lambda i: 1.0, range(1,len(titles),2))
                ))
    
    def get_sensitivity_pannel(self):
        if not self.is_resistance_unit():
            return {}
        titles = self.get_childnode_titles()
        if len(self.scores)==len(titles):
            return dict(zip(
                        map(lambda i: titles[i].split("~")[-1], range(1,len(titles),2)),
                        map(lambda i: self.scores[i], range(0,len(self.scores),2))
                    ))
        # Empty terminal unit with antibiotic resistance set by default
        return dict(zip(
                    map(lambda i: titles[i].split("~")[-1], range(1,len(titles),2)),
                    map(lambda i: 0.0, range(0,len(titles),2))
                ))
    
    def get_childnode_titles(self):
        return map(lambda s: s.strip(), self.child_node_titles.split(" | "))

    def format_child_node_titles(self):
        if not self.child_node_titles:
            return ""
        child_nodes = self.child_node_titles.split(" | ")
        if len(child_nodes) == 1:
            return child_nodes[0]
        return " | ".join([child_nodes[0]]+map(lambda i: child_nodes[i], range(1,len(child_nodes),2)))
    
    def normalize_score(self,i):
        norm = sum(map(lambda oME: oME.entry.get_power(i), filter(lambda record: record.entry.pair == int(math.ceil((i+1)/2.0)), self)))
        if norm:
            return float(self.scores[i])/norm
        return 0
    
    def split_entries(self):
        if not len(self):
            return [self.copy()]
        pair_num = max(map(lambda oME: oME.entry.pair, self))
        if pair_num == 1:
            return [self.copy()]
        subset = []
        child_nodes = self.child_node_titles.split(" | ")
        for i in range(pair_num):
            items = filter(lambda oME: oME.entry.pair == i+1, self)
            oMU = MessengerUnit(self.title," | ".join(child_nodes[0+2*i:2+2*i]),len(items))
            oMU.scores = self.scores[0+2*i:2+2*i]
            oMU.descriptions = [self.descriptions[i]]
            oMU.power_values = [self.power_values[i]]
            oMU.set(items)
            subset.append(oMU.copy())
        return subset
    
    def enroll(self):
        output = []
        for oME in self:
            output.append("%d\t%s [%d]\t%s [%s]\t%s/%s;" % (oME.entry.position,oME.entry.locus,
                    oME.entry.codon,oME.allele,oME.entry.alleles,
                    format_number(oME.entry.get_score(0,oME.allele)),
                    format_number(oME.entry.get_score(1,oME.allele))))
        output.append("\t"*5 + "/".join(map(lambda v: format_number(v), self.scores)))
        return output
    
    def get_sequence(self):
        seq = ""
        for entry in self:
            seq += entry.allele
        return seq
    
    def copy(self):
        oMU = MessengerUnit(self.title,self.child_node_titles,self.size)
        oMU.child_node_titles = self.child_node_titles
        oMU.mismatch_threshold = self.mismatch_threshold
        oMU.coverage_threshold = self.coverage_threshold
        oMU.foldchange_threshold = self.foldchange_threshold
        oMU.input_type = self.input_type
        oMU.e_threshold = self.e_threshold
        oMU.score_threshold = self.score_threshold
        oMU.success_cutoff = self.success_cutoff
        oMU.scores = []
        oMU.scores.extend(self.scores)
        oMU.descriptions = []
        oMU.descriptions.extend(self.descriptions)
        oMU.power_values = []
        oMU.power_values.extend(self.power_values)
        oMU.clade = self.clade
        oMU.set(map(lambda oME: oME.copy(),self))
        return oMU

class MessengerEntry:
    def __init__(self,title,oKeyEntry,allele):
        self.title = title
        self.entry = oKeyEntry
        self.allele = allele
        
    def get_maxScore(self):
        return self.entry.get_maxScore(self.allele)
                
    def copy(self):
        return MessengerEntry(self.title,self.entry.copy(),self.allele)
    
class Scores(collection.Collection):    # collection of ScoreUnits
    def __init__(self):
        collection.Collection.__init__(self)
        self.clade = ""
        self.winner = None
        
    def set_scores(self,scores):
        self.set(scores)
        self.set_winner()
        
    def set_winner(self):
        if not len(self):
            return
        self.sort()
        self.winner = self[0]
        if len(self) > 1:
            self.winner += self[1]
        self.clade = self.winner.title
        
    def sort(self):
        if len(self) < 2:
            return
        self.container.sort(self._sort)
        
    def _sort(self,oSU1,oSU2):
        return -cmp(oSU1.score,oSU2.score)
    
    def tostring(self,keys=[]):
        if self.winner == None:
            self.set_winner()
        return "Clade: " + oSU.tostring(key)

###############################################################################
class ScoreUnit:
    def __init__(self,title,score,sites={},sensitivity_pannel={},resistance_pannel={},repr=""):
        self.title = title
        self.score = score
        self.sites = sites
        self.sensitivity_pannel = sensitivity_pannel
        self.resistance_pannel = resistance_pannel
        if not repr:
            self.repr = "%s (%s)" % (self.title,format_number(self.score))
        
    def __repr__(self):
        return self.repr
    
    def __add__(self,other):
        oSU = self.copy()
        oSU.repr = "%s || %s" % (self.repr,other.repr)
        for key in other.resistance_pannel.keys():
            if key in oSU.resistance_pannel and oSU.resistance_pannel[key] > other.resistance_pannel[key]:
                continue
            oSU.resistance_pannel[key] = other.resistance_pannel[key]
            oSU.sensitivity_pannel[key] = other.sensitivity_pannel[key]
            oSU.sites[key] = other.sites[key]
        return oSU
    
    def get_resistance_pannel(self):
        pannel = {}
        pannel.update(self.resistance_pannel)
        return pannel
    
    def get_sensitivity_pannel(self):
        pannel = {}
        pannel.update(self.sensitivity_pannel)
        return pannel
    
    def estimate_resistance(self,ant):
        if ant not in self.resistance_pannel and ant not in self.sensitivity_pannel:
            return [0,0]
        q = (1.0 + math.log((1.0+self.resistance_pannel[ant])/(1.0+self.sensitivity_pannel[ant]),2))/2.0
        err = 0.5
        if ant in self.sites and self.sites[ant] > 1:
            err = q*(1.0-q)/math.sqrt(self.sites[ant]-1)
        corr = 2.0**(2.0*self.resistance_pannel[ant]-1.0)
        res = q*corr
        if res > 1.0:
            res = 1.0
        err *= 2.0*corr
        return [res,err]
    
    def get_predicted_resistance(self,ls):
        pannel = {"Sens":self.predict_sensitivity(False)}
        for ant in ls:
            if ant in self.resistance_pannel:
                if not self.sites[ant]:
                    pannel[ant] = ["ND",0.0]
                else:
                    pannel[ant] = self.predict_resistance(ant,False)
            else:
                pannel[ant] = ["DEFAULT",0.0]
        return pannel
    
    def predict_resistance(self,ant,flg_tostring=True):
        if ant.upper() == "SENS":
            return self.predict_sensitivity(flg_tostring)
        if ant in self.resistance_pannel and ant in self.sensitivity_pannel:
            res,err = self.estimate_resistance(ant)
            if flg_tostring:
                return "%s +/- %s" % (format_number(res,2),format_number(err,2))
            return [res,err]
        return ""
    
    def predict_sensitivity(self,flg_tostring=True):
        ant,v = self.top_resistance()
        if not ant and not v:
            if flg_tostring:
                return "1.0 +/- 0.0"
            return [1.0,0.0]
        res,err = self.estimate_resistance(ant)
        if flg_tostring:
            return "%s +/- %s" % (format_number(1.0-res,2),format_number(err,2))
        return [1.0-res,err]
    
    def top_resistance(self):
        topres = ["",0]
        if not self.resistance_pannel:
            return topres
        for ant in self.resistance_pannel:
            if self.resistance_pannel[ant] > topres[1]:
                topres = [ant,self.resistance_pannel[ant]]
        return topres
    
    def legend(self,ant):
        if ant.upper()=="SENS":
            compound,v = self.top_resistance()
            if compound:
                info = "%d sites" % self.sites[compound]
            else:
                info = "0 sites"
            return "\t%s\t%s" % (ant,self.predict_resistance(ant))
        sens = self.sensitivity_pannel[ant]
        res = self.resistance_pannel[ant]
        if ant in self.sites and self.sites[ant]:
            info = "%d sites" % self.sites[ant]
        elif sens:
            info = "default"
        else:
            info = "0 sites"
        return "\t%s\t%s [%s]\t%s/%s" % (ant,self.predict_resistance(ant),info,
            format_number(sens),format_number(res))
    
    def tostring(self,keys = []):
        s = self.repr
        available_keys = self.resistance_pannel.keys()
        if self.resistance_pannel:
            s += "\nAntibiotic sensitivity/resistance:\n"
            if not keys:
                keys = available_keys
            s += "\n".join(map(lambda key: self.legend(key),filter(lambda k: k in available_keys, keys)+["Sens"]))
        return s
    
    def copy(self):
        return ScoreUnit(self.title,self.score,self.sites,self.sensitivity_pannel,self.resistance_pannel,self.repr)

class parent:
    def __init__(self,rootfolder,basefolder,path,oTable,oReference,moltype="CDS"):
        self.basefolder = basefolder
        self.rootfolder = os.path.join(self.basefolder,rootfolder)
        self.path = os.path.join(self.basefolder,path)
        self.moltype = moltype
        if not os.path.exists(self.path):
            self.path = ""
        if self.moltype == "CDS":
            self.source_path = os.path.join(self.rootfolder,"bin","tmp","tmp.faa")
        else:
            self.source_path = os.path.join(self.rootfolder,"bin","tmp","tmp.fnn")
        self.oTable = oTable
        self.oReference = oReference
        self.dbname = ""
        self.refgenome_length = 4411532
        
    def clean(self):
        path = os.path.join(self.basefolder,"bin","tmp")
        if os.path.exists(path):
            for fname in os.listdir(path):
                os.remove(os.path.join(path,fname))
            
    def __call__(self):
        return
    
class Genbank(parent):
    def __init__(self,rootfolder,basefolder,path,oTable,oReference,moltype="CDS"): # moltype CDS | gene
        parent.__init__(self,rootfolder,basefolder,path,oTable,oReference,moltype)
        
    def __call__(self):
        self._prepare_data()
        if self.moltype == "CDS":
            oMessenger = self.oTable.identify("CDS",os.path.join(self.rootfolder,"bin"),self.dbname,self.source_path,self.oReference)
        else:
            oMessenger = self.oTable.identify("contigs",os.path.join(self.rootfolder,"bin"),self.dbname,self.source_path,self.oReference)
        self.clean()
        return oMessenger
        
    def _prepare_data(self):
        cns_fasta = oIO.get_CDS_from_GBK(self.path)
        if os.path.exists(self.source_path):
            os.remove(self.source_path)
        oIO.save("\n".join(cns_fasta),self.source_path)
        self.dbname = blast.formatdb(os.path.join(self.basefolder,"bin"),
            self.source_path,os.path.join(self.rootfolder,"bin","tmp","db"))
    
class Consensus(parent):
    def __init__(self,rootfolder,basefolder,path,oTable,oReference):
        parent.__init__(self,rootfolder,basefolder,path,oTable,oReference)
        
    def __call__(self):
        query_fasta,path = oIO.openFasta(self.path)
        if len(query_fasta.values()[0]) != REFGENOME_LENGTH:
            obj = Contigs(self.rootfolder,self.basefolder,os.path.basename(self.path),self.oTable,self.oReference)
            return obj()
        oMessenger = self.oTable.identify("consensus",self.rootfolder,self.oReference,"",query_fasta)
        self.clean()
        return oMessenger
    
class Contigs(Genbank):
    def __init__(self,rootfolder,basefolder,path,oTable,oReference):
        Genbank.__init__(self,rootfolder,basefolder,path,oTable,oReference,"gene")
        self.mlength = 50000
        
    def _prepare_data(self):
        def check_fasta_file(fasta_dic):
            fasta = {}
            for key in fasta_dic.keys():
                length = len(fasta_dic[key])
                if length > self.mlength:
                    k = length/self.mlength+1
                    counter = 1
                    start = 0
                    step = length/k
                    stop = start + step + 1000
                    while stop < length:
                        fasta["%s|%d" % (key,counter)] = fasta_dic[key][start:stop]
                        start += step
                        stop = start + step + 1000
                        counter += 1
                    fasta["%s|%d" % (key,counter)] = fasta_dic[key][start:length]
                else:
                    fasta[key] = fasta_dic[key]
            return fasta
            
        self.clean()
        dna_fasta,path = oIO.openFasta(self.path)
        fasta = check_fasta_file(dna_fasta)
        oIO.save("\n".join(map(lambda key: ">%s\n%s" % (key,fasta[key]),fasta.keys())),self.source_path)
        self.dbname = blast.formatdb(os.path.join(self.rootfolder,"bin"),
            self.source_path,os.path.join(self.rootfolder,"bin","tmp","db"),False)

class Protein(Genbank):
    def __init__(self,rootfolder,basefolder,path,oTable,oReference):
        Genbank.__init__(self,rootfolder,basefolder,path,oTable,oReference)
        
    def _prepare_data(self):
        def sort_seq(a,b):
            return cmp(int(a.split(" | ")[0]),int(b.split(" | ")[0]))
        self.clean()
        cns_fasta,path = oIO.openFasta(self.path)
        seqnames = cns_fasta.keys()
        seqnames.sort(sort_seq)
        
        oIO.save("\n".join(map(lambda key: ">%s\n%s" % (key,cns_fasta[key]),seqnames)),self.source_path)
        self.dbname = blast.formatdb(os.path.join(self.rootfolder,"bin"),
            self.source_path,os.path.join(self.rootfolder,"bin","tmp","db"))
    
class Gene(Genbank):
    def __init__(self,rootfolder,basefolder,path,oTable,oReference):
        Genbank.__init__(self,rootfolder,basefolder,path,oTable,oReference)
        
    def _prepare_data(self):
        def sort_seq(a,b):
            return cmp(int(a.split(" | ")[0]),int(b.split(" | ")[0]))
        self.clean()
        dna_fasta,path = oIO.openFasta(self.path)
        seqnames = dna_fasta.keys()
        seqnames.sort(sort_seq)
        dna_fasta = dict(zip(seqnames,oIO.translate(map(lambda key: dna_fasta[key],seqnames))))
        oIO.save("\n".join(map(lambda key: ">%s\n%s" % (key,dna_fasta[key]),seqnames)),self.source_path)
        self.dbname = blast.formatdb(os.path.join(self.rootfolder,"bin"),
            self.source_path,os.path.join(self.rootfolder,"bin","tmp","db"))
    
class VCF(Genbank):
    def __init__(self,rootfolder,basefolder,path,oTable,oReference):
        parent.__init__(self,rootfolder,basefolder,path,oTable,oReference)
        
    def __call__(self):
        dataset = self._format_vcf(oIO.open_text_file(self.path,True,True))
        oMessenger = self.oTable.identify("vcf",self.rootfolder,dataset,"",self._format_ref())
        self.clean()
        return oMessenger
    
    def _format_vcf(self,dataset):
        dataset = map(lambda line: line.split("\t"), filter(lambda item: item and item[0] != "#", dataset))
        dataset = map(lambda item: [int(item[1]),item[3].upper(),item[4].upper(),float(item[5])], dataset)
        return dataset
    
    def _format_ref(self):
        def parse_key(key):
            dataset = key.split(" | ")
            coords = dataset[-1][1:-1].split("..")
            return [dataset[1],int(coords[0]),int(coords[1])]
        return map(lambda key: parse_key(key)+[self.oReference[key]], self.oReference.keys())
    
###############################################################################
if __name__ == "__main__":
    options = {
               "-i":"input",    # input folder
               "-o":"output",   # output folder
               "-f":"Mycobacterium_bovis" ,        # parent folder
               "-g":"",          # generic name of the reference file
               "-b":"..",       # base folder
               "-r":"",       # root folder
            }
    options["-g"] = filter(lambda f: f[-4:].upper()==".FAA", os.listdir(os.path.join(options['-b'],options['-f'],"sources")))[0][:-4]
    execute(options)
