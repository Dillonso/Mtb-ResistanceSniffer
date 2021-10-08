import sys, os, string, math, subprocess, re
import seq_io, StringIO

##############################################################################################################
class container:
    def __init__(self):
        self.container = []
    
    def __len__(self):
        if not self.container:
            return 0
        return len(self.container)
    
    def __getitem__(self,key):
        if type(key)==type(0):
            if len(self) <= key:
                return
            return self.container[key]
        elif type(key)==type(""):
            for record in self.container:
                if record.title == key:
                    return record
        return 
    
    def __contains__(self,key):
        for record in self.container:
            if record.title == key:
                return True
        return False
    
    def __iter__(self):
        if not self.container:
            return iter([])
        records = []
        for record in self.container:
            records.append(record.copy())
        return iter(records)

    # Container items have to have the method copy()
    def get(self):
        container = []
        for record in self.container:
            container.append(record.copy())
        return container

    def format_title(self,title,length,space): # space - number of characters per line
        words = ("%s (%d bp)" % (title,length)).split(" ")
        title = [""]
        for word in words:
            title[-1] += word+" "
            if len(title[-1]) >= space:
                title.append("")
        return title
                    
##############################################################################################################
class BLAST(container):
    ###################################################
    ####  program in blastn,blastp,bl2seq,bl2seqp
    ####  path - executables
    ####  ref - path, sequence, blast database or fasta
    ####        file to be converted to blastdb
    ####  query - path or sequence
    ###################################################
    def __init__(self,program,seqtype="dna",path="",ref="",query="",source_path=""): # seqtype in "dna","protein","codon","genome"
        self.path = path
        self.program = program
        self.seqtype = seqtype
        self.success = True
        self.IO = seq_io.IO()
        self.query = query
        self.query_title = ""
        self.query_length = 0
        self.query_genemap = []
        self.sbjct = ref
        self.sbjct_title = ""
        self.sbjct_length = 0
        self.sbjct_genemap = []
        self.cline = ""
        container.__init__(self)
        
        # CONSTANTS
        self.tmpname1 = self.IO.random_filename("???_?????.tmp",self.IO.listdir(self.path,"tmp"))
        self.tmpname2 = self.IO.random_filename("???_?????.tmp",self.IO.listdir(self.path,"tmp"))
        self.tmpname3 = self.IO.random_filename("???_?????.tmp",self.IO.listdir(self.path,"tmp"))
        
        self._validate(source_path)
        if self.success:
            self._set_cline()
    
    def execute(self,flg_print=False):
        if not self.success or not self.cline:
            return
        output = self.process(self.cline)
        if flg_print:
            print self.cline
            print output
        if not output:
            self.IO.clean(self.path,[self.tmpname1,self.tmpname2,self.tmpname3])
            return False
        oParser = parser(self.program,self.seqtype,output)
        self.container = oParser()
        self.IO.clean(self.path,[self.tmpname1,self.tmpname2,self.tmpname3])
        return True
    
    def process(self,cline,flg_wait=False):
        try:
            process = subprocess.Popen(cline,stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                shell=(sys.platform!="win32"))
            if flg_wait:
                process.wait()
            process.stdin.close()
            output = process.stdout.read()
            process.stdout.close()
            return output
        except:
            return None
    
    def tostring(self,e_threshold=None):
        query_summary = header = output = ""
        sbjct_summary = {}
        for record in self:
            query_name = record.title
            header += "Query: " + query_name + "\n\n"
            output += "Query: " + query_name + "\n\n"
            for i in range(len(record.descriptions)):
                info = record.descriptions[i]
                if e_threshold != "None" and e_threshold != None and e_threshold < info.e:
                    continue
                alignment = record[i]
                header += "\tSubject: " + alignment.source + ("\t%f\t%f" % (info.score,info.e)) + "\n"
                if alignment.source not in sbjct_summary:
                    sbjct_summary[alignment.source] = []
                sbjct_summary[alignment.source].append(query_name)
                output += ("\tSubject: "+alignment.title+"\n\n")
                for hsp in alignment:
                    if hsp.positives == None:
                        output += ("\tScore = %f, E-value = %f,\n\tIdentities = %f\n\n" %
                            (hsp.score,hsp.expect,hsp.identities))
                    else:
                        output += ("\tScore = %f, E-value = %f,\n\tIdentities = %f, Positives = %f\n\n" %
                            (hsp.score,hsp.expect,hsp.identities,hsp.positives))
                    point = 0
                    q_gaps = s_gaps = 0
                    while point < hsp.alignment_length:
                        end = point+60
                        if end > hsp.alignment_length:
                            end = hsp.alignment_length
                        q_substr = hsp.query[point:end].upper()
                        s_substr = hsp.sbjct[point:end].upper()
                        qs = q_gaps
                        ss = s_gaps
                        q_gaps += q_substr.count("-")
                        s_gaps += s_substr.count("-")
                        output += "%g\t\t%s\t\t%g\n" % (hsp.query_start+point-qs,q_substr,hsp.query_start+end-1-q_gaps)
                        output += "%g\t\t%s\t\t%g\n\n" % (hsp.sbjct_start+point-ss,s_substr,hsp.sbjct_start+end-1-s_gaps)
                        point += 60
                output += "\n"
            header += "\n"
        header += "#"*60
        if sbjct_summary:
            query_summary += "Subject summary:\n"
            sbjct_summary = sbjct_summary.items()
            if len(sbjct_summary) > 1:
                sbjct_summary.sort(self.sortfn1)
                sbjct_summary.reverse()
            for item in sbjct_summary:
                query_summary += "\t".join(["",str(len(item[1])),item[0]])+"\n"
        return "\n\n".join([query_summary,header,output])
    
    def svg(self,X=25,Y=25,width=900,height=200,flg_finish=True):
        title_width=200
        title_height = 50
        title_start = width-title_width+5
        if not len(self):
            return ""
        max_length = max(self.query_length,self.sbjct_length)
        c = float(width-title_width)/max_length
        query_svg = self.svg_genemap(X,Y,self.query_title,self.query_length,self.query_genemap,
            title_width+c*self.query_length,title_height,title_width,title_start,False)
        span = height-2*title_height+20
        sbjct_svg = self.svg_genemap(X,Y+span,self.sbjct_title,self.sbjct_length,self.sbjct_genemap,
            title_width+c*self.sbjct_length,title_height,title_width,title_start,False)
        svg = [query_svg,sbjct_svg]
        if len(self) == 1 and self[0][0]:
            alignment_summary = self[0][0].summarize()
        else:
            alignment_summary = self.summarize()
        svg.append("<text x=\"%d\" y=\"%d\">Program - %s; Sequence type - %s; Summarized score = %d; Best expectation = %f</text>" % 
            (X,Y,self.program,self.seqtype,alignment_summary.score,alignment_summary.expect))
        zoom = 1.0
        for record in self:
            for sbjct in record:
                if self.seqtype == "genome":
                    zoom = 3.0
                    hsp_svg = sbjct.svg_alignment(c,record.shift,sbjct.shift,zoom,100,
                        X,Y+title_height,width-title_width,span-title_height,False)
                else:
                    hsp_svg = sbjct.svg_hsps(c,record.shift,sbjct.shift,zoom,50,
                        X,Y+title_height,width-title_width,span-title_height,False)
                if hsp_svg:
                    svg.append(hsp_svg)
        if flg_finish:
            svg.insert(0,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,Y+600))
            svg.append("</svg>")
        return "\n".join(svg)
    
    def svg_genemap(self,X=5,Y=25,genome_title="",genome_length=0,genemap=[],width=800,height=50,title_width=100,title_start=0,flg_finish=True):
        svg = []
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%f\" />" %
            (X,Y+height/2,X+width-title_width,Y+height/2,"red",3.0))
        # GENE MAP
        c = float(width-title_width)/genome_length
        for gene in genemap:
            try:
                lb,rb = map(lambda s: int(s),gene.split("-"))
            except:
                lb,rb = map(lambda s: int(s),gene.split(".."))
            if lb < 0:
                lb = 0
            if rb > genome_length:
                rb = genome_length
            if (genemap[gene]['remark'].find("hypothetical") > -1 or
                genemap[gene]['name'].find("hypothetical") > -1 or
                genemap[gene]['remark'].find("unknown") > -1 or
                genemap[gene]['name'].find("unknown") > -1):
                color = "grey"
            else:
                color = "green"
            bar_height = height/5
            shift = height/7
            if genemap[gene]['direction'] == "rev":
                shift = height-bar_height-height/7
            if not genemap[gene]['name']:
                title = genemap[gene]['remark']
            elif not genemap[gene]['remark']:
                title = genemap[gene]['name']
            else:
                title = genemap[gene]['name'] + "; " + genemap[gene]['remark']
            svg.append("<a title=\"%s\"><rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:%s;stroke:%s\" /></a>" %
                (title,X+c*lb,Y+shift,c*(rb-lb),bar_height,color,"grey"))
        # TITLE
        title = self.format_title(genome_title,genome_length,title_width/5)
        for subtitle in title:
            svg.append("<text x=\"%d\" y=\"%d\">%s</text>" % (X+title_start,Y+5,subtitle))
            Y += 15
        if flg_finish:
            svg.insert(0,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,height))
            svg.append("</svg>")
        return "\n".join(svg)
    
    def get_record(self,title):
        for record in self.container:
            if record.title == title:
                return record.copy()
        return None
        
    def get_top_alignment(self,title="locus"):
        for record in self:
            if record.title == title and record.container:
                hsps = []
                for hsp in record.container[0]:
                    hsps.append(hsp.copy())
                return [record.container[0].e,record.container[0].score,record.container[0].title,hsps]
        return None
    
    def get_close_to_point(self,point,e_threshold=0):
        hits = []
        for record in self.container:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect > e_threshold:
                        continue
                    hits.append([point-hsp.query_start,
                        point-hsp.query_end,
                        alignment.title,hsp.expect])
        if len(hits) > 1:
            hits.sort(self.sort_neighbours)
        return hits

    def _set_cline(self,blastdb=None):
        if self.seqtype in ("dna","codon"):
            algorithm = "blastn"
            flg_prot = "F"
            dbtype = "nucl"
        elif self.seqtype in ("protein","genome"):
            algorithm = "blastp"
            flg_prot = "T"
            dbtype = "prot"
        if sys.platform == "win32":
            if blastdb != None:
                self.cline = "%s -i %s -p %s -n %s" % (os.path.join(self.path,"formatdb.exe"),self.sbjct,flg_prot,blastdb)
            elif self.program == "bl2seq" and self.seqtype != "genome":
                default = "-G 6 -E 2 -q -2 -r 1 -e 1.0 -F F"
                self.cline = "%s -p %s -i %s -j %s %s" % (os.path.join(self.path,"bl2seq"),algorithm,self.query,self.sbjct,default)
            elif self.program == "bl2seq" and self.seqtype == "genome":
                default = "-G 6 -E 2 -q -2 -r 1 -e 10.0"
                self.cline = "%s -p %s -i %s -d %s %s" % (os.path.join(self.path,"blastall"),algorithm,self.query,self.sbjct,default)
            elif self.program == "blast":
                default = "-G 6 -E 2 -q -2 -r 1 -e 1.0"
                self.cline = "%s -p %s -i %s -d %s %s" % (os.path.join(self.path,"blastall"),algorithm,self.query,self.sbjct,default)
        elif sys.platform == "linux2":
            if blastdb != None:
                self.cline = "%s -in %s -dbtype %s -out %s" % (os.path.join(self.path,"makeblastdb"),self.sbjct,dbtype,blastdb)
            elif self.program == "blast":
                default = "-gapopen 6 -gapextend 2 -evalue 1.0"
                self.cline = "%s -query %s -db %s %s" % (os.path.join(self.path,algorithm),self.query,self.sbjct,default)
            elif self.program == "bl2seq" and self.seqtype == "genome":
                default = "-gapopen 6 -gapextend 2 -evalue 10.0 -dust no"
                self.cline = "%s -query %s -db %s %s" % (os.path.join(self.path,algorithm),self.query,self.sbjct,default)
            elif self.program == "bl2seq":
                default = "-gapopen 6 -gapextend 2 -evalue 1.0"
                self.cline = "%s -query %s -subject %s %s" % (os.path.join(self.path,algorithm),self.query,self.sbjct,default)

    def _validate(self,source_path=""):
        if self.program == "bl2seq":
            if (self.seqtype == "genome" and 
                    (not self.query or not self.sbjct or 
                    self.query[self.query.rfind("."):].upper() not in (".GBK",".GB") or
                    self.sbjct[self.sbjct.rfind("."):].upper() not in (".GBK",".GB"))):
                    return
            # 1) query is a StringIO object
            if type(self.query) == type(StringIO.StringIO()):
                self.query = self.query.read()
            # 2) query is a sequence string - ">seqname\nsequence or sequence"
            if self.query and not os.path.exists(self.query):
                self.query = self.query.replace("\r","")
                try:
                    seqname,sequence = self.query.split("\n")
                except:
                    seqname = ">Seq1"
                    sequence = self.query
                self.IO.save("%s\n%s" % (seqname,sequence),os.path.join(self.path,self.tmpname1))
                self.query = os.path.join(self.path,self.tmpname1)
                self.query_length = len(sequence)
                self.query_title = seqname
            # 3) query is a GBK file
            elif self.query and os.path.exists(self.query) and self.query[self.query.rfind("."):].upper() in (".GBK",".GB"):
                if self.seqtype == "dna":
                    self.query = self.IO.save(self.IO.openGBK(self.query,"FASTA"),os.path.join(self.path,self.tmpname1))
                elif self.seqtype == "genome":
                    self.query = self.IO.save_genes2fasta(self.query,os.path.join(self.path,self.tmpname1),True)
                    if not self.query:
                        self.success = False
                        return
                else:
                    self.success = False
                    return
                self.query_genemap = self.IO.getGeneMap()
                self.query_length = len(self.IO.getSequence())
                self.query_title = self.IO.getName()
            if not os.path.exists(self.query):
                self.success = False
                return
            
            # 1) sbjct is a StringIO object
            if type(self.sbjct) == type(StringIO.StringIO()):
                self.sbjct = self.sbjct.read()
            # 2) sbjct is a sequence string - ">seqname\nsequence or sequence"
            if self.sbjct and not os.path.exists(self.sbjct):
                self.sbjct = self.sbjct.replace("\r","")
                try:
                    seqname,sequence = self.sbjct.split("\n")
                except:
                    seqname = ">Seq2"
                    sequence = self.sbjct
                self.IO.save("%s\n%s" % (seqname,sequence),os.path.join(self.path,self.tmpname2))
                self.sbjct = os.path.join(self.path,self.tmpname2)
                self.sbjct_length = len(sequence)
                self.sbjct_title = seqname
            # 3) query is a GBK file
            elif self.sbjct and os.path.exists(self.sbjct) and self.sbjct[self.sbjct.rfind("."):].upper() in (".GBK",".GB"):
                if self.seqtype == "dna":
                    self.sbjct = self.IO.save(self.IO.openGBK(self.sbjct,"FASTA"),os.path.join(self.path,self.tmpname2))
                elif self.seqtype == "genome":
                    self.sbjct = self.IO.save_genes2fasta(self.sbjct,os.path.join(self.path,self.tmpname2),True)
                    self.sbjct = self.fasta2blast()
                else:
                    self.success = False
                    return
                self.sbjct_genemap = self.IO.getGeneMap()
                self.sbjct_length = len(self.IO.getSequence())
                self.sbjct_title = self.IO.getName()
            if not self.sbjct:
                self.success = False
                return
        else:
            # query
            if type(self.query) == type(StringIO.StringIO()):
                self.query = self.query.read()
            if self.query and not os.path.exists(self.query):
                self.IO.save(">locus\n"+self.query,os.path.join(self.path,self.tmpname2))
                self.query = os.path.join(self.path,self.tmpname2)
            elif self.query and self.query[self.query.rfind("."):].upper() in (".GBK",".GB"):
                outpath = os.path.join(self.path,self.tmpname2)
                self.query = os.path.join(source_path,os.path.basename(self.query))
                if self.seqtype == "dna":
                    flg_protein = False
                elif self.seqtype == "protein":
                    flg_protein = True
                self.query = self.IO.save_genes2fasta(self.query,outpath,flg_protein)
            if not os.path.exists(self.query):
                self.success = False
                return
            # subject
            if self.sbjct and self.sbjct[self.sbjct.rfind("."):].upper() in (".GBK",".GB"):
                outpath = os.path.join(self.path,self.tmpname1)
                self.sbjct = os.path.join(source_path,os.path.basename(self.sbjct))
                if self.program == "blastn":
                    flg_protein = False
                elif self.program == "blastp":
                    flg_protein = True
                self.sbjct = self.IO.save_genes2fasta(self.sbjct,outpath,flg_protein)
            elif self.sbjct and self.sbjct[self.sbjct.rfind("."):].upper() in (".FASTA",".FST",".FAS",".FNA",".FSA",".FAA",".FNN"):
                outpath = os.path.join(self.path,self.tmpname1)
                self.sbjct = os.path.join(source_path,os.path.basename(self.sbjct))
                self.sbjct = self.IO.copy(self.sbjct,outpath)
            if not self.sbjct:
                self.success = False
                return
            elif os.path.exists(self.sbjct):
                self.sbjct = self.fasta2blast()
                if not self.sbjct:
                    self.success = False
                    return
            else:
                pass

    def summarize(self):
        summarized_alignment = None
        for record in self:
            for alignment in record:
                summary = alignment.summarize()
                if not summarized_alignment:
                    summarized_alignment = summary.copy()
                else:
                    summarized_alignment += summary
        return summarized_alignment
    
    def fasta2blast(self):
        if not os.path.exists(self.sbjct):
            return
        blastdb = os.path.join(self.path,self.tmpname3)
        self._set_cline(blastdb)
        output = self.process(self.cline,True)
        if output == None:
            return
        return blastdb
    
    def sortfn1(self,a,b):
        return len(a[1]) - len(b[1])

    def sortfn2(self,a,b):
        return int(b[1]) - int(a[1])

    def sort_neighbours(self,a,b):
        return min(abs(a[0]),abs(a[1]))-min(abs(b[0]),abs(b[1]))
    
##############################################################################################################
class blast_record(container):
    def __init__(self,query_length,title=""):
        self.query_length = query_length
        self.title = title
        self.top_score = 0
        self.summerized_score = 0.0
        self.descriptions = []
        self.shift = 0
        container.__init__(self)
        
    def copy(self):
        record = blast_record(self.query_length,self.title)
        record.top_score = self.top_score
        record.summerized_score = self.summerized_score
        record.shift = self.shift
        for description in self.descriptions:
            record.add_description(description.title,description.score,description.e)
        for alignment in self.container:
            record.add_alignment(alignment.copy())
        return record
        
    def add_description(self,title,score,e):
        description = blast_description(title,score,e)
        self.descriptions.append(description)
    
    def add_alignment(self,alignment):
        self.container.append(alignment)
    
    def svg(self,X=5,Y=25,genemap=[],width=800,height=50,title_width=100,title_start=0,flg_finish=True):
        svg = []
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%f\" />" %
            (X,Y+height/2,X+width-title_width,Y+height/2,"red",3.0))
        # GENE MAP
        c = float(width-title_width)/self.query_length
        for gene in genemap:
            try:
                lb,rb = map(lambda s: int(s),gene.split("-"))
            except:
                lb,rb = map(lambda s: int(s),gene.split(".."))
            if lb < 0:
                lb = 0
            if rb > self.query_length:
                rb = self.query_length
            if (genemap[gene]['remark'].find("hypothetical") > -1 or
                genemap[gene]['name'].find("hypothetical") > -1 or
                genemap[gene]['remark'].find("unknown") > -1 or
                genemap[gene]['name'].find("unknown") > -1):
                color = "grey"
            else:
                color = "green"
            bar_height = height/5
            shift = height/7
            if genemap[gene]['direction'] == "rev":
                shift = height-bar_height-height/7
            if not genemap[gene]['name']:
                title = genemap[gene]['remark']
            elif not genemap[gene]['remark']:
                title = genemap[gene]['name']
            else:
                title = genemap[gene]['name'] + "; " + genemap[gene]['remark']
            svg.append("<a title=\"%s\"><rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:%s;stroke:%s\" /></a>" %
                (title,X+c*lb,Y+shift,c*(rb-lb),bar_height,color,"grey"))
        # TITLE
        title = self.format_title(self.title,self.query_length,title_width/5)
        for subtitle in title:
            svg.append("<text x=\"%d\" y=\"%d\">%s</text>" % (X+title_start,Y+5,subtitle))
            Y += 15
        if flg_finish:
            svg.insert(0,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,height))
            svg.append("</svg>")
        return "\n".join(svg)
    
    def sort(self):
        if len(self.container) > 1:
            self.container.sort(self._alignment_sort)

    def _alignment_sort(self,a,b):
        if a.e == b.e:
            return int(b.score-a.score)
        return int(a.e-b.e)

##############################################################################################################
class blast_description:
    def __init__(self,title,score,e):
        self.title = title
        self.score = score
        self.e = e
        
##############################################################################################################
class blast_alignment(container):
    def __init__(self,sbjct_length=0,title="",score=0,e=None,hsps=[]):
        self.sbjct_length = sbjct_length
        self.title = title
        self.source = self.parse_title()
        self.score = score
        self.e = e
        self.shift = 0
        container.__init__(self)
        self.container.extend(hsps)

    def copy(self):
        hsps = []
        for hsp in self.container:
            hsps.append(hsp.copy())
        alignment = blast_alignment(self.sbjct_length,self.title,self.score,self.e,hsps)
        alignment.shift = self.shift
        return alignment
    
    def __add__(self,other=None):
        new_alignment = self.copy()
        if not other:
            return new_alignment
        if new_alignment.title:
            new_alignment.title += "; "+other.title
        else:
            new_alignment.title = other.title
        if new_alignment.source:
            new_alignment.source += "; "+other.source
        else:
            new_alignment.source = other.source
        if new_alignment.score:
            new_alignment.score += other.score
        else:
            new_alignment.score = other.score
        if new_alignment.e != None:
            new_alignment.e = min(new_alignment.e,other.e)
        else:
            new_alignment.e = other.e
        new_alignment.shift = min(new_alignment.shift,other.shift)
        new_alignment.hsps.extend(other.hsps)
        return new_alignment

    def summarize(self):
        if not len(self.container):
            return None
        summarized_alignment = self.container[0].copy()
        for i in range(1,len(self.container)):
            summarized_alignment += self.container[i]
        return summarized_alignment
    
    def get_score(self):
        if not self.container:
            return None
        score = 0
        for hsp in self.container:
            score += hsp.score
        return score
    
    def get_expect(self):
        if not self.container:
            return None
        if len(self.container) == 1:
            return self.container[0].expect
        expect = self.container[0].expect
        for i in range(len(self.container),1):
            if self.container[i].expect < expect:
                expect = self.container[i].expect
        return expect
    
    def svg(self,X=5,Y=25,genemap=[],width=800,height=50,title_width=100,title_start=0,flg_finish=True):
        svg = []
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%f\" />" %
            (X,Y+height/2,X+width-title_width,Y+height/2,"red",3.0))
        # GENE MAP
        c = float(width-title_width)/self.sbjct_length
        if not genemap:
            genemap = {}
        for gene in genemap:
            try:
                lb,rb = map(lambda s: int(s),gene.split("-"))
            except:
                lb,rb = map(lambda s: int(s),gene.split(".."))
            if lb < 0:
                lb = 0
            if rb > self.sbjct_length:
                rb = self.sbjct_length
            if (genemap[gene]['remark'].find("hypothetical") > -1 or
                genemap[gene]['name'].find("hypothetical") > -1 or
                genemap[gene]['remark'].find("unknown") > -1 or
                genemap[gene]['name'].find("unknown") > -1):
                color = "grey"
            else:
                color = "green"
            bar_height = height/5
            shift = height/7
            if genemap[gene]['direction'] == "rev":
                shift = height-bar_height-height/7
            if not genemap[gene]['name']:
                title = genemap[gene]['remark']
            elif not genemap[gene]['remark']:
                title = genemap[gene]['name']
            else:
                title = genemap[gene]['name'] + "; " + genemap[gene]['remark']
            svg.append("<a title=\"%s\"><rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:%s;stroke:%s\" /></a>" %
                (title,X+c*lb,Y+shift,c*(rb-lb),bar_height,color,"grey"))
        # TITLE
        title = self.format_title(self.title,self.sbjct_length,title_width/5)
        for subtitle in title:
            svg.append("<text x=\"%d\" y=\"%d\">%s</text>" % (X+title_start,Y+5,subtitle))
            Y += 15
        if flg_finish:
            svg.insert(0,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,height))
            svg.append("</svg>")
        return "\n".join(svg)
    
    def svg_hsps(self,c,query_shift,sbjct_shift,zoom=1.0,score_threshold=50,
                    X=5,Y=25,width=800,height=150,flg_finish=True):
        svg = []
        for hsp in self:
            if hsp.score < score_threshold:
                continue
            if not hsp.strand or hsp.strand=="Plus/Plus":
                svg.append("<path d=\"M%f,%dL%f,%dL%f,%dL%f,%dZ\" style=\"fill:%s;stroke:%s;opacity:%f\" />" %
                    (X+c*(query_shift+hsp.query_start*zoom),Y,X+c*(query_shift+hsp.query_end*zoom),Y,
                    X+c*(sbjct_shift+hsp.sbjct_end*zoom),Y+height,X+c*(sbjct_shift+hsp.sbjct_start*zoom),Y+height,
                    "blue","grey",0.5))
            elif hsp.strand=="Plus/Minus":
                svg.append("<path d=\"M%f,%dL%f,%dL%f,%dL%f,%dZ\" style=\"fill:%s;stroke:%s;opacity:%f\" />" %
                    (X+c*(query_shift+hsp.query_start*zoom),Y,X+c*(query_shift+hsp.query_end*zoom),Y,
                    X+c*(sbjct_shift+hsp.sbjct_start*zoom),Y+height,X+c*(sbjct_shift+hsp.sbjct_end*zoom),Y+height,
                    "blue","grey",0.5))
        if flg_finish:
            svg.insert(0,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,height))
            svg.append("</svg>")
        return "\n".join(svg)
    
    def svg_alignment(self,c,query_shift,sbjct_shift,zoom=1.0,score_threshold=100,
                        X=5,Y=25,width=800,height=150,flg_finish=True):
        svg = []
        query_end = sbjct_end = score = 0
        query_start = sbjct_start = self.sbjct_length
        for hsp in self:
            if hsp.query_end > query_end:
                query_end = hsp.query_end
            if hsp.sbjct_end > sbjct_end:
                sbjct_end = hsp.sbjct_end
            if hsp.query_start < query_start:
                query_start = hsp.query_start
            if hsp.sbjct_start < sbjct_start:
                sbjct_start = hsp.sbjct_start
            score += hsp.score
        if score < score_threshold:
            return ""
        svg.append("<a title=\"Score = %d\"><path d=\"M%f,%dL%f,%dL%f,%dL%f,%dZ\" style=\"fill:%s;stroke:%s;opacity:%f\" /></a>" %
            (score,X+c*(query_shift+query_start*zoom),Y,X+c*(query_shift+query_end*zoom),Y,
            X+c*(sbjct_shift+sbjct_end*zoom),Y+height,X+c*(sbjct_shift+sbjct_start*zoom),Y+height,
            "blue","none",score/1000.0))
        if flg_finish:
            svg.insert(0,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,height))
            svg.append("</svg>")
        return "\n".join(svg)
    
    def parse_title(self):
        if not self.title:
            return ""
        p = self.title.find("; from ")
        if p == -1:
            return self.title
        source = self.title.replace(", complete sequence.","")
        return source[p+7:]

##############################################################################################################
class blast_hsp:
    def __init__(self,score,e,aln_length,identities,positives,gaps,strand,qlb,slb,qrb,srb,query,sbjct,hits):
        self.score = score
        self.expect = e
        self.alignment_length = int(aln_length)
        self.identities = identities
        self.positives = positives
        self.gaps = gaps
        self.strand = strand
        self.query_start = qlb
        self.sbjct_start = slb
        self.query_end = qrb
        self.sbjct_end = srb
        self.query = query
        self.sbjct = sbjct
        self.match = hits
    
    def __str__(self):
        return "Query [%d..%d]; Sbjct [%d..%d]" % (self.query_start,self.query_end,self.sbjct_start,self.sbjct_end)
    
    def copy(self):
        hsp = blast_hsp(self.score,
                            self.expect,
                            self.alignment_length,
                            self.identities,
                            self.positives,
                            self.gaps,
                            self.strand,
                            self.query_start,
                            self.sbjct_start,
                            self.query_end,
                            self.sbjct_end,
                            self.query,
                            self.sbjct,
                            self.match)
        return hsp

    def __add__(self,other=None):
        new_hsp = self.copy()
        if not other:
            return new_hsp
        new_hsp.score += other.score
        new_hsp.expect = min(new_hsp.expect,other.expect)
        if other.identities:
            new_hsp.identities += other.identities
        if other.positives:
            new_hsp.positives += other.positives
        if other.gaps:
            new_hsp.gaps += other.gaps
        new_hsp.strand = ""
        new_hsp.query_start = min(new_hsp.query_start,other.query_start)
        new_hsp.query_stop = max(new_hsp.query_end,other.query_end)
        new_hsp.sbjct_start = min(new_hsp.sbjct_start,other.sbjct_start)
        new_hsp.sbjct_stop = max(new_hsp.sbjct_end,other.sbjct_end)
        new_hsp.query = ""
        new_hsp.sbjct = ""
        new_hsp.match = ""
        return new_hsp

##############################################################################################################
class parser:
    def __init__(self,program,seqtype,raw_text):
        if sys.platform == "win32":
            self.oParser = win_parser(program,seqtype,raw_text)
        elif sys.platform == "linux2":
            self.oParser = linux_parser(program,seqtype,raw_text)
        
    def __call__(self):
        return self.oParser._parse()
     
##############################################################################################################
class sys_parser:
    def __init__(self,program,seqtype,raw_text):
        self.program = program
        self.seqtype = seqtype
        self.raw_text = raw_text
        self.dataset = []
        
    def _parse(self):
        self.raw_text = self.raw_text.replace("\r","")
        '''
        #### TEMP
        IO = seq_io.IO()
        IO.save(self.raw_text,"wintmp.out","a")
        '''
        if self.program == "bl2seq" and self.seqtype != "genome":
            return self._parse_bl2seq()
        if self.program == "bl2seq" and self.seqtype == "genome":
            return self._parse_genome()
        q = self.raw_text.find("Query= ")
        if q == -1:
            return
        while q != -1:
            t = self.raw_text.find("Query= ",q+1)
            title = self.raw_text[q+7:self.raw_text.find("\n",q)]
            query_length = int(self.raw_text[self.raw_text.find("         (",q)+10:self.raw_text.find(" letters)\n",q)].replace(",",""))
            record = blast_record(query_length,title)
            a = self.raw_text.find("Sequences producing significant alignments",q,t)
            if a == -1:
                q = self.raw_text.find("Query= ",q+1)
                continue
            else:
                a = self.raw_text.find("\n\n",a)+2
            b = self.raw_text.find(">",a)-2
            headers = self.raw_text[a:b].split("\n")
            for hit in headers:
                if not hit:
                    continue
                e,score,name = self._parse_spsa(hit)
                if e == None:
                    continue
                record.add_description(name,score,e)
                alignment = self._parse_alignment(name,self.raw_text[q:t])
                if not alignment:
                    alignment = blast_alignment()
                record.add_alignment(alignment)
            self.dataset.append(record)
            q = self.raw_text.find("Query= ",q+1)
        return self.dataset
            
    def _parse_bl2seq(self):
        q = self.raw_text.find("Query= ")
        title = self.raw_text[q+7:self.raw_text.find("\n",q)]
        query_length = int(self.raw_text[self.raw_text.find("         (",q)+10:self.raw_text.find(" letters)\n",q)].replace(",",""))
        record = blast_record(query_length,title)
        a = self.raw_text.find(">")
        name = self.raw_text[a+1:self.raw_text.find("\n",a)]
        alignment = self._parse_alignment(name,self.raw_text)
        if alignment:
            record.add_description(name,alignment.get_score(),alignment.get_expect())
            record.add_alignment(alignment)
        self.dataset.append(record)
        return self.dataset
    
    def _parse_genome(self):
        q = self.raw_text.find("Query= ")
        if q == -1:
            return
        while q != -1:
            t = self.raw_text.find("Query= ",q+1)
            g = self.raw_text.find("\n         (",q)
            title = self.raw_text[q+7:g].replace("\n"," ")
            query_length = int(self.raw_text[g+11:self.raw_text.find(" letters)\n",q)].replace(",",""))
            record = blast_record(query_length,title)
            record.shift = self._get_shift(title)
            a = self.raw_text.find("Sequences producing significant alignments",q,t)
            if a == -1:
                q = self.raw_text.find("Query= ",q+1)
                continue
            else:
                a = self.raw_text.find("\n\n",a)+2
            b = self.raw_text.find(">",a)-2
            headers = self.raw_text[a:b].split("\n")
            for hit in headers:
                if not hit:
                    continue
                e,score,name = self._parse_spsa(hit)
                if e == None:
                    continue
                record.add_description(name,score,e)
                alignment = self._parse_alignment(name,self.raw_text[q:t])
                if not alignment:
                    alignment = blast_alignment()
                alignment.shift = self._get_shift(alignment.title)
                record.add_alignment(alignment)
            self.dataset.append(record)
            q = self.raw_text.find("Query= ",q+1)
        return self.dataset
            
    def _parse_spsa(self,hit):
        name = hit[:73]
        while name[-1] == " ":
            name = name[:-1]
        try:
            e = self._format_e(hit[73:])
            score = int(name[name.rfind(" ")+1:])
        except:
            return None,None,""
        name = name[:name.rfind(" ")]
        while name[-1] == " ":
            name = name[:-1]
        return e,score,name
            
    def _parse_alignment(self,name,input):
        if len(name) > 20:
            name = name[:21]
        p = input.find(">"+name)
        if p == -1:
            return
        l = input.find("          Length = ",p)
        title = input[p+1:l]
        sbjct_length = int(input[l+19:input.find("\n",l)].replace(",",""))
        for symbol in ["\n","\r","\t","  "]:
            title = title.replace(symbol,"")
        hsps = []
        d = input.find(">",p+1)
        if d == -1:
            d = input.find("  Database:",p+1)
        block = input[p:d]
        b = block.find(" Score = ")
        score = e = None
        while b > -1:
            q = block.find("Query:",b)
            score,e,aln_length,identities,positives,gaps,strand = self._parse_hsp_header(block[b:q])
            qlb,slb,qrb,srb,query,sbjct,hits = self._parse_hsp_body(block[q:block.find(" Score = ",b+1)])
            hsps.append(blast_hsp(score,e,aln_length,identities,positives,gaps,strand,qlb,slb,qrb,srb,query,sbjct,hits))
            b = block.find(" Score = ",b+1)
        if not score:
            return
        alignment = blast_alignment(sbjct_length,title,score,e,hsps)
        return alignment
    
    def _parse_hsp_header(self,block):
        s = block.find(" Score = ")
        m = block.find(",   Method:")
        d = block.find(" Identities = ")
        p = block.find(" Positives = ")
        g = block.find("Gaps = ",d)
        t = block.find("Strand")
        score = float(block[s+9:block.find(" bits ")])
        if m > -1:
            e = self._format_e(block[block.find("Expect = ")+9:m])
        else:
            try:
                e = self._format_e(block[block.find("Expect = ")+9:block.find("\n",s)])
            except:
                ex = block.find("Expect = ")+9
                e = self._format_e(block[ex:block.find(",",ex)])
        aln_length = int(block[block.find("/",d)+1:block.find(" (",d)])
        identities = int(block[d+13:block.find("/",d)])
        positives = None
        if p > -1:
            positives = int(block[p+12:block.find("/",p)])
        gaps = 0
        if g > -1:
            gaps = int(block[g+7:block.find("/",g)])
        strand = ""
        if t > -1:
            strand = block[block.find("Plus",t):block.find("\n",t)]
            strand = strand.replace(" ","")
        return score,e,aln_length,identities,positives,gaps,strand
        
    def _parse_hsp_body(self,block):
        block = block.replace("\r","")
        lines = filter(lambda l: l, block.split("\n"))
        query = sbjct = hits = ""
        i = j = 0
        while i < len(lines):
            if i == 0:
                qlb = int(lines[i][6:lines[i].find(" ",8)])
                slb = int(lines[i+2][6:lines[i+2].find(" ",8)])
            if len(lines[i]) < 12 or lines[i][:5] != "Query":
                i += 1
                continue
            indend = lines[i][:14].rfind(" ")+1
            dedend = lines[i].find(" ",indend+1)
            query += lines[i][indend:dedend]
            hits += lines[i+1][indend:dedend]
            sbjct += lines[i+2][indend:dedend]
            j = i
            i += 3
        qrb = int(lines[j][lines[j].rfind(" "):])
        srb = int(lines[j+2][lines[j+2].rfind(" "):])
        return qlb,slb,qrb,srb,query,sbjct,hits
    '''
    def _parse_hsp_body(self,block):
        lines = block.split("\n")
        query = sbjct = hits = ""
        i = j = 0
        while i < len(lines):
            if i == 0:
                qlb = int(lines[i][6:lines[i].find(" ",8)])
                slb = int(lines[i+2][6:lines[i+2].find(" ",8)])
            if len(lines[i]) < 12 or lines[i][:5] != "Query":
                i += 1
                continue
            indend = lines[i][:14].rfind(" ")+1
            dedend = lines[i].find(" ",13)
            query += lines[i][indend:dedend]
            hits += lines[i+1][indend:dedend]
            sbjct += lines[i+2][indend:dedend]
            j = i
            i += 3
        qrb = int(lines[j][lines[j].rfind(" "):])
        srb = int(lines[j+2][lines[j+2].rfind(" "):])
        return qlb,slb,qrb,srb,query,sbjct,hits
    '''
    def _get_shift(self,title):
        a = title.find("[")
        b = title.find("..",a)
        if a > -1 and b > -1 and title.find("]",b) > -1:
            try:
                return int(title[a+1:b])
            except:
                return 0
        return 0
        
    def _format_e(self,e):  # e as string;
        if e.find(".") > -1:
            e = float(e)
        else:
            d = e.find("e-")
            try:
                x = int(e[:d])
            except:
                x = 1
            y = int(e[d+2:])
            e = x*10.0**(-y)
        return e

##############################################################################################################
class win_parser(sys_parser):
    def __init__(self,program,seqtype,raw_text):
        sys_parser.__init__(self,program,seqtype,raw_text)

##############################################################################################################
class linux_parser(sys_parser):
    def __init__(self,program,seqtype,raw_text):
        sys_parser.__init__(self,program,seqtype,raw_text)

    def _parse(self):
        self.raw_text = self.raw_text.replace("\r","")
        #### TEMP
        IO = seq_io.IO()
        IO.save(self.raw_text,"lintmp.out")
        if self.program in ("bl2seq","bl2seqp") and self.seqtype != "genome":
            return self._parse_bl2seq()
        q = self.raw_text.find("Query= ")
        if q == -1:
            return
        while q != -1:
            t = self.raw_text.find("Query= ",q+1)
            l = self.raw_text.find("Length=",q+1)+7
            title = self.raw_text[q+7:self.raw_text.find("\n",q)]
            query_length = int(self.raw_text[l:self.raw_text.find("\n",l)].replace(",",""))
            record = blast_record(query_length,title)
            a = self.raw_text.find("Sequences producing significant alignments",q,t)
            if a == -1:
                q = self.raw_text.find("Query= ",q+1)
                continue
            else:
                a = self.raw_text.find("\n\n",a)+2
            b = self.raw_text.find(">",a)-2
            headers = self.raw_text[a:b].split("\n")
            for hit in headers:
                if not hit:
                    continue
                result = self._parse_spsa(hit)
                e,score,name = result
                if e == None:
                    continue
                record.add_description(name,score,e)
                alignment = self._parse_alignment(name,self.raw_text[q:t])
                if not alignment:
                    alignment = blast_alignment()
                record.add_alignment(alignment)
            self.dataset.append(record)
            q = self.raw_text.find("Query= ",q+1)
        return self.dataset
            
    def _parse_bl2seq(self):
        q = self.raw_text.find("Query= ")
        l = self.raw_text.find("Length=",q+1)+7
        title = self.raw_text[q+7:self.raw_text.find("\n",q)]
        query_length = int(self.raw_text[l:self.raw_text.find("\n",l)].replace(",",""))
        record = blast_record(query_length,title)
        alignment = self._parse_bl2seq_alignment()
        if alignment:
            record.add_description(alignment.title,alignment.get_score(),alignment.get_expect())
            record.add_alignment(alignment)
            self.dataset.append(record)
        return self.dataset
    
    def _parse_spsa(self,hit):
        name = hit[2:78]
        while name[-1] == " ":
            name = name[:-1]
        try:
            e = self._format_e(hit[78:])
            score = float(name[name.rfind(" ")+1:])
        except:
            return None,None,""
        name = name[:name.rfind(" ")]
        while name[-1] == " ":
            name = name[:-1]
        return e,score,name
            
    def _parse_alignment(self,name,input):
        if len(name) > 20:
            name = name[:21]
        p = input.find("> "+name)
        if p == -1:
            return
        l = input.find("Length=",p)+7
        title = input[p+2:input.find("\n",p)]
        sbjct_length = int(input[l:input.find("\n",l)].replace(",",""))
        for symbol in ["\n","\r","\t","  "]:
            title = title.replace(symbol,"")
        hsps = []
        d = input.find(">",p+1)
        if d == -1:
            d = input.find("  Database:",p+1)
        block = input[p:d]
        b = block.find(" Score = ")
        while b > -1:
            q = block.find("Query",b)
            score,e,aln_length,identities,positives,gaps,strand = self._parse_hsp_header(block[b:q])
            qlb,slb,qrb,srb,query,sbjct,hits = self._parse_hsp_body(block[q:block.find(" Score = ",b+1)])
            hsps.append(blast_hsp(score,e,aln_length,identities,positives,gaps,strand,qlb,slb,qrb,srb,query,sbjct,hits))
            b = block.find(" Score = ",b+1)
        alignment = blast_alignment(sbjct_length,title,score,e,hsps)
        return alignment
            
    def _parse_bl2seq_alignment(self):
        p = self.raw_text.find("Subject=")+8
        title = self.raw_text[p:self.raw_text.find("\n",p)]
        l = self.raw_text.find("Length=",p)+7
        sbjct_length = int(self.raw_text[l:self.raw_text.find("\n",l)].replace(",",""))
        for symbol in ["\n","\r","\t","  "]:
            title = title.replace(symbol,"")
        hsps = []
        block = self.raw_text[p:]
        b = block.find(" Score = ")
        while b > -1:
            q = block.find("Query",b)
            score,e,aln_length,identities,positives,gaps,strand = self._parse_hsp_header(block[b:q])
            b = block.find(" Score = ",b+1)
            qlb,slb,qrb,srb,query,sbjct,hits = self._parse_hsp_body(block[q:b])
            hsps.append(blast_hsp(score,e,aln_length,identities,positives,gaps,strand,qlb,slb,qrb,srb,query,sbjct,hits))
        try:
            alignment = blast_alignment(sbjct_length,title,score,e,hsps)
        except:
            alignment = None
        return alignment
    '''
    def _parse_hsp_body(self,block):
        block = block.replace("\r","")
        lines = filter(lambda l: l, block.split("\n"))
        query = sbjct = hits = ""
        i = j = 0
        while i < len(lines):
            if i == 0:
                qlb = int(lines[i][6:lines[i].find(" ",8)])
                slb = int(lines[i+2][6:lines[i+2].find(" ",8)])
            if len(lines[i]) < 12 or lines[i][:5] != "Query":
                i += 1
                continue
            indend = lines[i][:14].rfind(" ")+1
            dedend = lines[i].find(" ",indend+1)
            query += lines[i][indend:dedend]
            hits += lines[i+1][indend:dedend]
            sbjct += lines[i+2][indend:dedend]
            j = i
            i += 3
        qrb = int(lines[j][lines[j].rfind(" "):])
        srb = int(lines[j+2][lines[j+2].rfind(" "):])
        return qlb,slb,qrb,srb,query,sbjct,hits
    '''
def formatdb(bin,infile,dbname,flg_prot=True,flg_wait=False):
    if not flg_prot:
        flg_prot = "F"
        dbtype = "nucl"
    else:
        flg_prot = "T"
        dbtype = "prot"
    if sys.platform == "win32":
        cline = "%s -i %s -p %s -n %s" % (os.path.join(bin,"formatdb.exe"),infile,flg_prot,dbname)
    elif sys.platform == "linux2":
        cline = "%s -in %s -dbtype %s -out %s" % (os.path.join(bin,"makeblastdb"),infile,dbtype,dbname)
    
    try:
        process = subprocess.Popen(cline,stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            shell=(sys.platform!="win32"))
        if flg_wait:
            process.wait()
        process.stdin.close()
        output = process.stdout.read()
        process.stdout.close()
        return dbname
    except:
        return infile
