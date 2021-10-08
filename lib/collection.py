import string

##############################################################################################################
class container:
    def __init__(self,name="",id=0):
        self.id = id
        self.title = name
        self.container = []
        self.flg_unique_entries = False
        self.level = 0
        self.trigger = None
        
    #### TEMP
    def check_trigger(self,mytitle=''):
        responce = []
        mytitle = "%s|%s" % (mytitle,self.title)
        for item in self:
            try:
                responce += item.check_trigger(mytitle)
            except:
                pass
        if self.trigger:
            print "collect:23",self.trigger
            try: 
                return [mytitle]
            except:
                return ['inidentified']
        return []
    
    def reset_trigger(self):
        self.trigger = None
        for item in self:
            try:
                item.reset_trigger()
            except:
                break
    
    #### Methods
    def __len__(self):
        if not self.container:
            return 0
        return len(self.container)
    
    def _append(self,item,index=None):
        if item == None:
            return False
        if self.flg_unique_entries and hasattr(item,'title') and item.title in self:
            return False
        if index == None or index < 0 or index >= len(self):
            self.container.append(item)
        else:
            self.container.insert[index,item]
        return True
    
    def _get_key(self,key):
        if type(key)==type(0):
            if len(self) <= key:
                return
            return key
        elif type(key)==type(""):
            try:
                return self.container.index(key)
            except:
                pass
            try:
                return map(lambda obj: self._get_title(obj) == key,self.container).index(True)
            except:
                return
    
    def _get_title(self,obj):
        if hasattr(obj, 'title'):
            return obj.title
        return ""

    def __getitem__(self,key):
        key = self._get_key(key)
        if key != None:
            return self.container[key]
    
    def __contains__(self,key):
        if self._get_key(key) != None:
            return True
        return False
    
    def __iter__(self):
        if not self.container:
            return iter([])
        records = []
        for record in self.container:
            records.append(record)
        return iter(records)
    
    def __setitem__(self,key,value):
        key = self._get_key(key)
        if key != None:
            self.container[key] = value
    
    def __delitem__(self,key):
        key = self._get_key(key)
        if key != None:
            del self.container[key]
    
    def __repr__(self):
        return ",".join(map(lambda item:str(item),self.get_titles()))
    
    # Sort methods:
    def _sort_by_id(self,a,b):
        return cmp(a.id,b.id)
    
    def _sort_by_title(self,a,b):
        return cmp(a.title,b.title)
    
    def append(self,item):
        return self._append(item)
    
    def extend(self,values):
        map(self._append,values)

    def insert(self,item,index=0):
        return self._append(item,index)
    
    def replace(self,ind,item):
        self[ind] = item
    
    def pop(self):
        if len(self):
            try:
                obj = self[-1].copy()
            except:
                return
            del self[-1]
            return obj
        return

    def has(self,key):
        try:
            return key in map(lambda item: self._get_title(item), self.container)
        except:
            return
        
    def index(self,title):
        if type(title) == type(0):
            index = title
            if index < -1 or index >= len(self):
                return
            return index
        titles = self.get_titles()
        if title in titles:
            return titles.index(title)
        return
    
    def set(self,data):
        self.container = []
        self.container.extend(data)

    # Container items have to have the method copy()
    def get(self,key=None,flg_make_copy=True):
        if key != None and type(key) == type(True):
            flg_make_copy = key
            key = None
        if key != None:
            item = self.__getitem__(key)
            if item != None and flg_make_copy:
                return item.copy()
            elif item != None and not flg_make_copy:
                return item
            else:
                return
        container = []
        for record in self.container:
            if flg_make_copy:
                try:
                    container.append(record.copy())
                except:
                    container.append(record)
            else:
                container.append(record)
        return container
    
    def get_titles(self):
        titles = map(lambda item: self._get_title(item),self.container)
        if len(titles) > 1:
            titles.sort()
        return titles
    
    def get_dictionary(self):
        if len(self):
            keys = self.get_titles()
            return dict(zip(keys, map(self.get,keys)))
        return {}

    def clear(self):
        self.id = 0
        self.title = ""
        self.container = []
        self.level = 0
        
    def copy(self):
        return self

##############################################################################
class stack(container):
    def __init__(self,name="",id=0):
        container.__init__(self,name,id)
        
    def pop(self):
        if not len(self.container):
            return
        obj = self.container[-1]
        del self.container[-1]
        return obj
    
    def push(self,obj):
        self.append(obj)
    
    def remove_all(self):
        self.container = []

##############################################################################
class Collection(container):
    def __init__(self,name="",id=0):
        container.__init__(self,name,id)
        self.level = 0

##############################################################################
class RootNode(container):
    def __init__(self,name="",id=0):
        container.__init__(self,name,id)
        self.level = 0

##############################################################################
class Node(container):
    def __init__(self,name="",id=0):
        container.__init__(self,name,id)
        self.level = 1

##############################################################################
class UndoRedo:
    def __init__(self):
        self.undo_stack = stack("Undo")
        self.redo_stack = stack("Redo")
    
    def push(self,obj):
        self.undo_stack.push(obj)
        self.redo_stack.remove_all()
        
    def undo(self):
        obj = self.undo_stack.pop()
        if obj:
            self.redo_stack.push(obj)
        return obj
    
    def redo(self):
        obj = self.redo_stack.pop()
        if obj:
            self.undo_stack.push(obj)
        return obj

###############################################################################

if __name__ == "__main__":
    oMain = Main("gene_repositary")
    oMain.execute()
