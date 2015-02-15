#!/usr/bin/python

"""\
html_writer.py - Construct HTML pages

"""

import os
import sys
import util
from subprocess import Popen

browser_bin = None
for file in ['C:\Program Files\Safari\Safari.exe', \
             '/Applications/Safari.app/Contents/MacOS/Safari',
             '/usr/bin/firefox']:
    if (os.path.isfile(file)):
        browser_bin = file
        print >> sys.stderr, "Using Browser: %s" % browser_bin
        break

class HtmlWriter:
    def __init__(self, filename, force_path_creation=True):
        self.filename = filename
        self.filepath = os.path.dirname(filename)
        if (not os.path.exists(self.filepath)):
            if (force_path_creation):
                util._mkdir(self.filepath)
            else:
                raise Exception("cannot write to HTML file %s since the directory doesn't exist" % filename)
        
        self.file = open(self.filename, "w")
        self.write("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" ")
        self.write("\"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">")
        self.write("<html>\n<body>\n")
        
        self.abspath = os.path.abspath(self.filename)
        return
    
    def get_abspath(self):
        return self.abspath
    
    def write(self, str):
        if (self.file == None):
            raise Exception("cannot write to this HTML since it is already closed")
        self.file.write(str)
    
    def link(self, url, text):
        self.write('<a href="%s">%s</a>' % (url, text))        
    
    def flush(self):
        if (self.file != None):
            self.file.flush()
    
    def write_svg(self, scene, relative_path):
        scene.embed_in_html(self.file, self.filepath, relative_path)
    
    def embed_svg(self, url, width, height):
        self.write('<embed src="%s" width="%d" height="%d" type="image/svg+xml">' % (url, width, height))
    
    def branch(self, relative_path, link_text=None):
        """Branches the HTML file by creating a new HTML and adding a link to it with the desired text
        """
        if (link_text == None):
            link_text = relative_path
            
        self.write("<a href=\"" + relative_path + ".html\">" + link_text + "</a>")
        self.flush()
        return HtmlWriter(os.path.join(self.filepath, relative_path + ".html"))
    
    def close(self):
        self.write("</body>\n</html>\n")
        self.file.flush()
        self.file.close()
        self.file = None
        
    def display(self):
        self.close()
        if (browser_bin != None):
            p = Popen([browser_bin, os.path.abspath(self.filename)])
            os.waitpid(p.pid, 0)
        else:
            print >> sys.stderr, "Cannot find a browser, won't display HTMLs"    

def test():
    html_write = HtmlWriter("../results/test.html")
    html_write.write("hello world")
    html_write.display()
    return

if __name__ == '__main__': test()
