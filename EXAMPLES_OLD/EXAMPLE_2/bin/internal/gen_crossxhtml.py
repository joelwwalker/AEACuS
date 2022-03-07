################################################################################
#
# Copyright (c) 2011 The MadGraph Development team and Contributors
#
# This file is a part of the MadGraph 5 project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph license which should accompany this 
# distribution.
#
# For more information, please visit: http://madgraph.phys.ucl.ac.be
#
################################################################################
""" Create gen_crossxhtml """


import os
import math
import re
import pickle
import re
try:
    import internal.files as files
    import internal.save_load_object as save_load_object
except:
    import madgraph.iolibs.files as files
    import madgraph.iolibs.save_load_object as save_load_object

pjoin = os.path.join
exists = os.path.exists



crossxhtml_template = """
<HTML> 
<HEAD> 
    %(refresh)s 
    <META HTTP-EQUIV="EXPIRES" CONTENT="20" > 
    <TITLE>Online Event Generation</TITLE>
    <link rel=stylesheet href="./HTML/mgstyle.css" type="text/css">
</HEAD>
<BODY>
<script type="text/javascript">
function UrlExists(url) {
  var http = new XMLHttpRequest();
  http.open('HEAD', url, false);
  try{
     http.send()
     }
  catch(err){
   return 1==2;
  }
  return http.status!=404;
}
function check_link(url,alt, id){
    var obj = document.getElementById(id);
    if ( ! UrlExists(url)){
       if ( ! UrlExists(alt)){
         obj.href = url;
         return 1==1;
        }
       obj.href = alt;
       return 1 == 2;
    }
    obj.href = url;
    return 1==1;
}
</script>    
    <H2 align=center> Results in the %(model)s for %(process)s </H2> 
    <HR>
    %(status)s
    <br>
    <br>
    <H2 align="center"> Available Results </H2>
        <TABLE BORDER=2 align="center">  
            <TR align="center">
                <TH>Run</TH> 
                <TH>Collider</TH> 
                <TH> Banner </TH>
                <TH> %(numerical_title)s </TH> 
                <TH> Events  </TH>
                <TH> Data </TH>  
                <TH>Output</TH>
                <TH>Action</TH> 
            </TR>      
            %(old_run)s
        </TABLE>
    <H3 align=center><A HREF="./index.html"> Main Page </A></H3>
</BODY> 
</HTML> 
"""

status_template = """
<H2 ALIGN=CENTER> Currently Running </H2>
<TABLE BORDER=2 ALIGN=CENTER>
    <TR ALIGN=CENTER>
        <TH nowrap ROWSPAN=2 font color="#0000FF"> Run Name </TH>
        <TH nowrap ROWSPAN=2 font color="#0000FF"> Tag Name </TH>
        <TH nowrap ROWSPAN=2 font color="#0000FF"> Cards </TH>   
        <TH nowrap ROWSPAN=2 font color="#0000FF"> Results </TH> 
        <TH nowrap ROWSPAN=1 COLSPAN=3 font color="#0000FF"> Status/Jobs </TH>
    </TR>
        <TR> 
            <TH>   Queued </TH>
            <TH>  Running </TH>
            <TH> Done  </TH>
        </TR>
    <TR ALIGN=CENTER> 
        <TD nowrap ROWSPAN=2> %(run_name)s </TD>
        <TD nowrap ROWSPAN=2> %(tag_name)s </TD>
        <TD nowrap ROWSPAN=2> <a href="./Cards/param_card.dat">param_card</a><BR>
                    <a href="./Cards/run_card.dat">run_card</a><BR>
                    %(plot_card)s
                    %(pythia_card)s
                    %(pgs_card)s
                    %(delphes_card)s
        </TD>
        <TD nowrap ROWSPAN=2> %(results)s </TD> 
        %(status)s
 </TR>
 <TR></TR>
   %(stop_form)s
 </TABLE>
"""

class AllResults(dict):
    """Store the results for all the run of a given directory"""
    
    web = False 
    
    def __init__(self, model, process, path):
        
        dict.__init__(self)
        self.order = []
        self.lastrun = None
        self.process = ', '.join(process)
        if len(self.process) > 60:
            pos = self.process[50:].find(',')
            if pos != -1:
                self.process = self.process[:50+pos] + ', ...'
        self.path = path
        self.model = model
        self.status = ''
        self.unit = 'pb'
        self.current = None
    
    def def_current(self, run, tag=None):
        """define the name of the current run
            The first argument can be a OneTagResults
        """

        if isinstance(run, OneTagResults):
            self.current = run
            self.lastrun = run['run_name']
            return
        
        assert run in self or run == None
        self.lastrun = run
        if run:
            if not tag:
                self.current = self[run][-1]
            else:
                self.current = self[run][tag]
        else:
            self.current = None
    
    def delete_run(self, run_name, tag=None):
        """delete a run from the database"""

        assert run_name in self

        if not tag :
            if self.current and self.current['run_name'] == run_name:
                self.def_current(None)                    
            del self[run_name]
            self.order.remove(run_name)
            if self.lastrun == run_name:
                self.lastrun = None
        else:
            assert tag in [a['tag'] for a in self[run_name]]
            RUN = self[run_name]
            if len(RUN) == 1:
                self.delete_run(run_name)
                return
            RUN.remove(tag)

        #update the html
        self.output()
    
    def def_web_mode(self, web):
        """define if we are in web mode or not """
        if web is True:
            try:
                web = os.environ['SERVER_NAME']
            except:
                web = 'my_computer'
        self['web'] = web
        self.web = web
        
    def add_run(self, name, run_card, current=True):
        """ Adding a run to this directory"""
        
        tag = run_card['run_tag']
        if name in self.order:
            #self.order.remove(name) # Reorder the run to put this one at the end 
            if  tag in self[name].tags:
                if self[name].return_tag(tag).parton and len(self[name]) > 1:
                    #move the parton information before the removr
                    self[name].return_tag(self[name][1]['tag']).parton = \
                                               self[name].return_tag(tag).parton
                if len(self[name]) > 1:        
                    self[name].remove(tag) # Remove previous tag if define 
                    self[name].add(OneTagResults(name, run_card, self.path))
            else:
                #add the new tag run    
                self[name].add(OneTagResults(name, run_card, self.path))
            new = self[name] 
        else:
            new = RunResults(name, run_card, self.process, self.path)
            self[name] = new  
            self.order.append(name)
        
        if current:
            self.def_current(name)        
        if new.info['unit'] == 'GeV':
            self.unit = 'GeV'
            
    def update(self, status, level, makehtml=True, error=False):
        """update the current run status"""
        if self.current:
            self.current.update_status(level)
        self.status = status
        if self.current and self.current.debug  and self.status and not error:
            self.current.debug = None

        if makehtml:
            self.output()

    def resetall(self, main_path=None):
        """check the output status of all run
           main_path redefines the path associated to the run (allowing to move 
           the directory)
        """
        
        self.path = main_path
        
        for key,run in self.items():
            if key == 'web':
                continue
            for i,subrun in enumerate(run):
                self.def_current(subrun)
                self.clean()
                self.current.event_path = pjoin(main_path,'Events') 
                if i==0:
                    self.current.update_status()
                else:
                    self.current.update_status(nolevel='parton')
        self.output()
                    
    def clean(self, levels = ['all'], run=None, tag=None):
        """clean the run for the levels"""

        if not run and not self.current:
            return
        to_clean = self.current
        if run and not tag:
            for tagrun in self[run]:
                self.clean(levels, run, tagrun['tag'])
            return

        if run:
            to_clean = self[run].return_tag(tag)
        else:
            run = to_clean['run_name']
        
        if 'all' in levels:
            levels = ['parton', 'pythia', 'pgs', 'delphes', 'channel']
        
        if 'parton' in levels:
            to_clean.parton = []
        if 'pythia' in levels:
            to_clean.pythia = []
        if 'pgs' in levels:
            to_clean.pgs = []
        if 'delphes' in levels:
            to_clean.delphes = []
        
        
    def save(self):
        """Save the results of this directory in a pickle file"""
        filename = pjoin(self.path, 'HTML', 'results.pkl')
        save_load_object.save_to_file(filename, self)

    def add_detail(self, name, value, run=None, tag=None):
        """ add information to current run (cross/error/event)"""
        assert name in ['cross', 'error', 'nb_event', 'cross_pythia',
                        'nb_event_pythia','error_pythia']

        if not run and not self.current:
            return

        if not run:
            run = self.current
        else:
            run = self[run].return_tag(tag)
            
        if name == 'cross_pythia':
            run['cross_pythia'] = float(value)
        elif name == 'nb_event':
            run[name] = int(value)
        elif name == 'nb_event_pythia':
            run[name] = int(value)
        else:    
            run[name] = float(value)    
    
    def output(self):
        """ write the output file """
        
        # 1) Create the text for the status directory        
        if self.status and self.current:
            if isinstance(self.status, str):
                status = '<td ROWSPAN=2 colspan=4>%s</td>' %  self.status
            else:
                s = self.status
                status ='''<td> %s </td> <td> %s </td> <td> %s </td>
                </tr><tr><td colspan=3><center> %s </center></td>''' % (s[0],s[1], s[2], s[3])
                
            
            status_dict = {'status': status,
                            'cross': self.current['cross'],
                            'error': self.current['error'],
                            'run_name': self.current['run_name'],
                            'tag_name': self.current['tag'],
                            'unit': self[self.current['run_name']].info['unit']}

            if exists(pjoin(self.path, 'HTML',self.current['run_name'], 
                        'results.html')):
                status_dict['results'] = """<A HREF="./HTML/%(run_name)s/results.html">%(cross).4g <font face=symbol>&#177;</font> %(error).4g (%(unit)s)</A>""" % status_dict
            else:
                status_dict['results'] = "No results yet"
            if exists(pjoin(self.path, 'Cards', 'plot_card.dat')):
                status_dict['plot_card'] = """ <a href="./Cards/plot_card.dat">plot_card</a><BR>"""
            else:
                status_dict['plot_card'] = ""
            if exists(pjoin(self.path, 'Cards', 'pythia_card.dat')):
                status_dict['pythia_card'] = """ <a href="./Cards/pythia_card.dat">pythia_card</a><BR>"""
            else:
                status_dict['pythia_card'] = ""
            if exists(pjoin(self.path, 'Cards', 'pgs_card.dat')):
                status_dict['pgs_card'] = """ <a href="./Cards/pgs_card.dat">pgs_card</a><BR>"""
            else:
                status_dict['pgs_card'] = ""
            if exists(pjoin(self.path, 'Cards', 'delphes_card.dat')):
                status_dict['delphes_card'] = """ <a href="./Cards/delphes_card.dat">delphes_card</a><BR>"""
            else:
                status_dict['delphes_card'] = ""
                
            if self.web:
                status_dict['stop_form'] = """
                 <TR ALIGN=CENTER><TD COLSPAN=7 text-align=center>
<FORM ACTION="http://%(web)s/cgi-bin/RunProcess/handle_runs-pl"  ENCTYPE="multipart/form-data" METHOD="POST">
<INPUT TYPE=HIDDEN NAME=directory VALUE="%(me_dir)s">
<INPUT TYPE=HIDDEN NAME=whattodo VALUE="stop_job">
<INPUT TYPE=SUBMIT VALUE="Stop Current Job">
</FORM></TD></TR>""" % {'me_dir': self.path, 'web': self.web}
            else:
                status_dict['stop_form'] = ""
            
            
            status = status_template % status_dict
            refresh = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"10\">"
        else:
            status =''
            refresh = ''
        
        
        # See if we need to incorporate the button for submission
        if os.path.exists(pjoin(self.path, 'RunWeb')):       
            running  = True
        else:
            running = False
        
        # 2) Create the text for the old run:
        old_run = ''
        for key in self.order:
            old_run += self[key].get_html(self.path, web=self.web, running=running)
        
        text_dict = {'process': self.process,
                     'model': self.model,
                     'status': status,
                     'old_run': old_run,
                     'refresh': refresh,
                     'numerical_title': self.unit == 'pb' and 'Cross section (pb)'\
                                                          or 'Width (GeV)'}
        
        text = crossxhtml_template % text_dict
        open(pjoin(self.path,'crossx.html'),'w').write(text)
        
       

class RunResults(list):
    """The list of all OneTagResults"""        

    def __init__(self, run_name, run_card, process, path):
        """initialize the object"""
        
        self.info = {'run_name': run_name,'me_dir':path}
        self.tags = [run_card['run_tag']]
        # Set the collider information
        data = process.split('>',1)[0].split()
        if len(data) == 2:
            name1,name2 = data
            if run_card['lpp1'] == '-1':
                name1 = ' p~'
            elif run_card['lpp1']  == '1':
                name1 = ' p'   
            elif run_card['lpp1'] == '2':
                name1 = ' a'
            if run_card['lpp2'] == '-1':
                name2 = 'p~'
            elif run_card['lpp2']  == '1':
                name2 = ' p' 
            elif run_card['lpp2'] == '2':
                name2 = ' a'                
            self.info['collider'] = '''%s %s <br> %s x %s  GeV''' % \
                    (name1, name2, run_card['ebeam1'], run_card['ebeam2'])
            self.info['unit'] = 'pb'                       
        else:
            self.info['collider'] = 'decay'
            self.info['unit'] = 'GeV'
        
        self.append(OneTagResults(run_name, run_card, path))
        
    
    def get_html(self, output_path, **opt):
        """WRITE HTML OUTPUT"""

        try:
            self.web = opt['web']
            self.info['web'] = self.web
        except:
            self.web = False

        # check if more than one parton output
        parton = [r for r in self if r.parton]
        # clean wrong previous run link
        if len(parton)>1:
            for p in parton[:-1]:
                p.parton = []

        dico = self.info
        dico['run_span'] = sum([tag.get_nb_line() for tag in self], 1) -1
        dico['tag_data'] = '\n'.join([tag.get_html(self) for tag in self])
        text = """
        <tr>
        <td rowspan=%(run_span)s>%(run_name)s</td> 
        <td rowspan=%(run_span)s><center> %(collider)s </center></td>
        %(tag_data)s
        </tr>
        """ % dico

        if self.web:
            
            text = text % self.info


        return text
    
    def return_tag(self, name):
        
        for data in self:
            if data['tag'] == name:
                return data
        
        if name is None:
            # return last entry
            return self[-1]
        
        raise Exception, '%s is not a valid tag' % name
    
    def is_empty(self):
        """Check if this run contains smtg else than html information"""

        if not self:
            return True
        if len(self) > 1:
            return False
        
        data = self[0]
        if data.parton or data.pythia or data.pgs or data.delphes:
            return False
        else:
            return True
        
    def add(self, obj):
        """ """
        
        assert isinstance(obj, OneTagResults)
        tag = obj['tag']
        assert tag not in self.tags
        self.tags.append(tag)
        self.append(obj)
        
    def get_last_pythia(self):
        for i in range(1, len(self)+1):
            if self[-i].pythia:
                return self[-i]['tag']

    def get_current_info(self):
        
        output = {}
        current = self[-1]
        # Check that cross/nb_event/error are define
        if current.pythia and not current['nb_event'] and len(self) > 1:
            output['nb_event'] = self[-2]['nb_event']
            output['cross'] = self[-2]['cross']
            output['error'] = self[-2]['error']
        elif (current.pgs or current.delphes) and not current['nb_event'] and len(self) > 1:
            if self[-2]['cross_pythia'] and self[-2]['nb_event_pythia']:
                output['cross'] = self[-2]['cross_pythia']
                output['nb_event'] = self[-2]['nb_event_pythia']
                output['error'] = self[-2]['error_pythia']
            else:
                output['nb_event'] = self[-2]['nb_event']
                output['cross'] = self[-2]['cross']
                output['error'] = self[-2]['error']
        elif current['cross']:
            return current
        elif len(self) > 1:
            output['nb_event'] = self[-2]['nb_event']
            output['cross'] = self[-2]['cross']
            output['error'] = self[-2]['error']
        else:
            output['nb_event'] = 0
            output['cross'] = 0
            output['error'] = 1e-99             
        return output
        
        
    def remove(self, tag):
        
        assert tag in self.tags
        
        obj = [o for o in self if o['tag']==tag][0]
        self.tags.remove(tag)
        list.remove(self, obj)
    
    
        
class OneTagResults(dict):
    """ Store the results of a specific run """
    
    def __init__(self, run_name, run_card, path):
        """initialize the object"""
        
        # define at run_result
        self['run_name'] = run_name
        self['tag'] = run_card['run_tag']
        self.event_path = pjoin(path,'Events')
        self.me_dir = path
        self.debug = None
        
        # Default value
        self['nb_event'] = 0
        self['cross'] = 0
        self['cross_pythia'] = ''
        self['nb_event_pythia'] = 0
        self['error'] = 0
        self.parton = [] 
        self.pythia = []
        self.pgs = []
        self.delphes = []
        # data 
        self.status = ''
        
    
    
    
    def update_status(self, level='all', nolevel=[]):
        """update the status of the current run """

        exists = os.path.exists
        run = self['run_name']
        tag =self['tag']
        
        path = pjoin(self.event_path, run)
        html_path = pjoin(self.event_path, os.pardir, 'HTML', run)
        
        # Check if the output of the last status exists
        if level in ['gridpack','all']:
            if 'gridpack' not in self.parton and \
                    exists(pjoin(path,os.pardir ,os.pardir,"%s_gridpack.tar.gz" % run)):
                self.parton.append('gridpack')
        
        if level in ['parton','all'] and 'parton' not in nolevel:
            
            if 'lhe' not in self.parton and \
                        (exists(pjoin(path,"unweighted_events.lhe.gz")) or
                         exists(pjoin(path,"unweighted_events.lhe"))):
                self.parton.append('lhe')
        
            if 'root' not in self.parton and \
                          exists(pjoin(path,"unweighted_events.root")):
                self.parton.append('root')
            
            if 'plot' not in self.parton and \
                                      exists(pjoin(html_path,"plots_parton.html")):
                self.parton.append('plot')

            if 'param_card' not in self.parton and \
                                    exists(pjoin(path, "param_card.dat")):
                self.parton.append('param_card')
                
        if level in ['pythia', 'all']:
            
            if 'plot' not in self.pythia and \
                          exists(pjoin(html_path,"plots_pythia_%s.html" % tag)):
                self.pythia.append('plot')
            
            if 'lhe' not in self.pythia and \
                            (exists(pjoin(path,"%s_pythia_events.lhe.gz" % tag)) or
                             exists(pjoin(path,"%s_pythia_events.lhe" % tag))):
                self.pythia.append('lhe')


            if 'hep' not in self.pythia and \
                            (exists(pjoin(path,"%s_pythia_events.hep.gz" % tag)) or
                             exists(pjoin(path,"%s_pythia_events.hep" % tag))):
                self.pythia.append('hep')
            
            if 'root' not in self.pythia and \
                              exists(pjoin(path,"%s_pythia_events.root" % tag)):
                self.pythia.append('root')
                
            if 'lheroot' not in self.pythia and \
                          exists(pjoin(path,"%s_pythia_lhe_events.root" % tag)):
                self.pythia.append('lheroot')
            
            if 'log' not in self.pythia and \
                          exists(pjoin(path,"%s_pythia.log" % tag)):
                self.pythia.append('log')     

        if level in ['pgs', 'all']:
            
            if 'plot' not in self.pgs and \
                         exists(pjoin(html_path,"plots_pgs_%s.html" % tag)):
                self.pgs.append('plot')
            
            if 'lhco' not in self.pgs and \
                              (exists(pjoin(path,"%s_pgs_events.lhco.gz" % tag)) or
                              exists(pjoin(path,"%s_pgs_events.lhco." % tag))):
                self.pgs.append('lhco')
                
            if 'root' not in self.pgs and \
                                 exists(pjoin(path,"%s_pgs_events.root" % tag)):
                self.pgs.append('root')
            
            if 'log' not in self.pgs and \
                          exists(pjoin(path,"%s_pgs.log" % tag)):
                self.pgs.append('log') 
    
        if level in ['delphes', 'all']:
            
            if 'plot' not in self.delphes and \
                              exists(pjoin(html_path,"plots_delphes_%s.html" % tag)):
                self.delphes.append('plot')
            
            if 'lhco' not in self.delphes and \
                 (exists(pjoin(path,"%s_delphes_events.lhco.gz" % tag)) or
                 exists(pjoin(path,"%s_delphes_events.lhco" % tag))):
                self.delphes.append('lhco')
                
            if 'root' not in self.delphes and \
                             exists(pjoin(path,"%s_delphes_events.root" % tag)):
                self.delphes.append('root')     
            
            if 'log' not in self.delphes and \
                          exists(pjoin(path,"%s_delphes.log" % tag)):
                self.delphes.append('log') 
    
    def special_link(self, link, level, name):
        
        id = '%s_%s_%s_%s' % (self['run_name'],self['tag'], level, name)
        
        return " <a  id='%(id)s' href='%(link)s.gz' onClick=\"check_link('%(link)s.gz','%(link)s','%(id)s')\">%(name)s</a>" \
              % {'link': link, 'id': id, 'name':name}
    
    def double_link(self, link1, link2, name, id):
        
         return " <a  id='%(id)s' href='%(link1)s' onClick=\"check_link('%(link1)s','%(link2)s','%(id)s')\">%(name)s</a>" \
              % {'link1': link1, 'link2':link2, 'id': id, 'name':name}       
        
    def get_links(self, level):
        """ Get the links for a given level"""
        
        out = ''
        if level == 'parton':
            if 'gridpack' in self.parton:
                out += self.special_link("./%(run_name)s_gridpack.tar",
                                                         'gridpack', 'gridpack')
            
            if 'lhe' in self.parton:
                link = './Events/%(run_name)s/unweighted_events.lhe'
                level = 'parton'
                name = 'LHE'
                out += self.special_link(link, level, name) 
            if 'root' in self.parton:
                out += ' <a href="./Events/%(run_name)s/unweighted_events.root">rootfile</a>'
            if 'plot' in self.parton:
                out += ' <a href="./HTML/%(run_name)s/plots_parton.html">plots</a>'
            if 'param_card' in self.parton:
                out += ' <a href="./Events/%(run_name)s/param_card.dat">param_card</a>'

            return out % self
        
        if level == 'pythia':          
            if 'log' in self.pythia:
                out += """ <a href="./Events/%(run_name)s/%(tag)s_pythia.log">LOG</a>"""
            if 'hep' in self.pythia:
                link = './Events/%(run_name)s/%(tag)s_pythia_events.hep'
                level = 'pythia'
                name = 'STDHEP'
                out += self.special_link(link, level, name)                 
            if 'lhe' in self.pythia:
                link = './Events/%(run_name)s/%(tag)s_pythia_events.lhe'
                level = 'pythia'
                name = 'LHE'                
                out += self.special_link(link, level, name) 
            if 'root' in self.pythia:
                out += """ <a href="./Events/%(run_name)s/%(tag)s_pythia_events.root">rootfile (LHE)</a>"""
            if 'lheroot' in self.pythia:
                out += """ <a href="./Events/%(run_name)s/%(tag)s_pythia_lhe_events.root">rootfile (LHE)</a>"""
            if 'plot' in self.pythia:
                out += ' <a href="./HTML/%(run_name)s/plots_pythia_%(tag)s.html">plots</a>'
            return out % self


        if level == 'pgs':
            if 'log' in self.pgs:
                out += """ <a href="./Events/%(run_name)s/%(tag)s_pgs.log">LOG</a>"""
            if 'lhco' in self.pgs:
                link = './Events/%(run_name)s/%(tag)s_pgs_events.lhco'
                level = 'pgs'
                name = 'LHCO'                
                out += self.special_link(link, level, name)  
            if 'root' in self.pgs:
                out += """ <a href="./Events/%(run_name)s/%(tag)s_pgs_events.root">rootfile</a>"""    
            if 'plot' in self.pgs:
                out += """ <a href="./HTML/%(run_name)s/plots_pgs_%(tag)s.html">plots</a>"""
            return out % self
        
        if level == 'delphes':
            
            if 'log' in self.delphes:
                out += """ <a href="./Events/%(run_name)s/%(tag)s_delphes.log">LOG</a>"""
            if 'lhco' in self.delphes:
                link = './Events/%(run_name)s/%(tag)s_delphes_events.lhco'
                level = 'delphes'
                name = 'LHCO'                
                out += self.special_link(link, level, name)
            if 'root' in self.delphes:
                out += """ <a href="./Events/%(run_name)s/%(tag)s_delphes_events.root">rootfile</a>"""    
            if 'plot' in self.delphes:
                out += """ <a href="./HTML/%(run_name)s/plots_delphes_%(tag)s.html">plots</a>"""            
            return out % self
                
    
    def get_nb_line(self):
        
        nb_line = 0
        for i in [self.parton, self.pythia, self.pgs, self.delphes]:
            if len(i):
                nb_line += 1
        return max([nb_line,1])
    
    
    def get_html(self, runresults):
        """create the html output linked to the this tag
           RunResults is given in case of cross-section need to be taken
           from a previous run
        """
        
        
        tag_template = """
        <td rowspan=%(tag_span)s> <a href="./Events/%(run)s/%(run)s_%(tag)s_banner.txt">%(tag)s</a>%(debug)s</td>
        %(subruns)s"""
        
        # Compute the text for eachsubpart
        
        sub_part_template_parton = """
        <td rowspan=%(cross_span)s><center><a href="./HTML/%(run)s/results.html"> %(cross).4g <font face=symbol>&#177;</font> %(err).2g </a></center></td>
        <td rowspan=%(cross_span)s><center> %(nb_event)s<center></td><td> %(type)s </td>
        <td> %(links)s</td>
        <td> %(action)s</td>
        </tr>"""
        
        sub_part_template_pgs = """
        <td> %(type)s </td>
        <td> %(links)s</td>
        <td> %(action)s</td> 
        </tr>"""        
        
        # Compute the HTMl output for subpart
        nb_line = self.get_nb_line()
        # Check that cross/nb_event/error are define
        if self.pythia and not self['nb_event']:
            try:
                self['nb_event'] = runresults[-2]['nb_event']
                self['cross'] = runresults[-2]['cross']
                self['error'] = runresults[-2]['error']
            except:
                pass
                
        elif (self.pgs or self.delphes) and not self['nb_event']:
            if runresults[-2]['cross_pythia'] and runresults[-2]['cross']:
                self['cross'] = runresults[-2]['cross_pythia']
                self['error'] = runresults[-2]['error_pythia']
                self['nb_event'] = runresults[-2]['nb_event_pythia']                           
            else:
                self['nb_event'] = runresults[-2]['nb_event']
                self['cross'] = runresults[-2]['cross']
                self['error'] = runresults[-2]['error']

        
        first = None
        subresults_html = ''
        for type in ['parton', 'pythia', 'pgs', 'delphes']:
            data = getattr(self, type)
            if not data:
                continue
            
            local_dico = {'type': type, 'run': self['run_name']}

            if not first:
                template = sub_part_template_parton
                first = type
                if type=='parton' and self['cross_pythia']:
                    local_dico['cross_span'] = 1
                    local_dico['cross'] = self['cross']
                    local_dico['err'] = self['error']
                    local_dico['nb_event'] = self['nb_event']
                elif self['cross_pythia']:
                    if self.parton:
                        local_dico['cross_span'] = nb_line -1
                    else:
                        local_dico['cross_span'] = nb_line
                    if self['nb_event_pythia']:
                        local_dico['nb_event'] = self['nb_event_pythia']
                    else:
                        local_dico['nb_event'] = 0
                    local_dico['cross'] = self['cross_pythia']
                    local_dico['err'] = self['error_pythia']
                else:
                    local_dico['cross_span'] = nb_line
                    local_dico['cross'] = self['cross']
                    local_dico['err'] = self['error']
                    local_dico['nb_event'] = self['nb_event']
                    
            elif type == 'pythia' and self['cross_pythia']:
                template = sub_part_template_parton
                if self.parton:           
                    local_dico['cross_span'] = nb_line - 1
                    if self['nb_event_pythia']:
                        local_dico['nb_event'] = self['nb_event_pythia']
                    else:
                        local_dico['nb_event'] = 0
                else:
                    local_dico['cross_span'] = nb_line
                    local_dico['nb_event'] = self['nb_event']
                local_dico['cross'] = self['cross_pythia']
                local_dico['err'] = self['error_pythia']
            else:
               template = sub_part_template_pgs             
            
            # Fill the links
            local_dico['links'] = self.get_links(type)

            # Fill the actions
            if type == 'parton':
                if runresults.web:
                    local_dico['action'] = """
<FORM ACTION="http://%(web)s/cgi-bin/RunProcess/handle_runs-pl"  ENCTYPE="multipart/form-data" METHOD="POST">
<INPUT TYPE=HIDDEN NAME=directory VALUE="%(me_dir)s">
<INPUT TYPE=HIDDEN NAME=whattodo VALUE="remove_level">
<INPUT TYPE=HIDDEN NAME=level VALUE="all">
<INPUT TYPE=HIDDEN NAME=tag VALUE=\"""" + self['tag'] + """\">
<INPUT TYPE=HIDDEN NAME=run VALUE="%(run_name)s">
<INPUT TYPE=SUBMIT VALUE="Remove run">
</FORM>
                    
<FORM ACTION="http://%(web)s/cgi-bin/RunProcess/handle_runs-pl"  ENCTYPE="multipart/form-data" METHOD="POST">
<INPUT TYPE=HIDDEN NAME=directory VALUE="%(me_dir)s">
<INPUT TYPE=HIDDEN NAME=whattodo VALUE="pythia">
<INPUT TYPE=HIDDEN NAME=run VALUE="%(run_name)s">
<INPUT TYPE=SUBMIT VALUE="Run Pythia">
</FORM>"""
                else:
                    local_dico['action'] = self.command_suggestion_html('remove %s parton --tag=%s' \
                                                                       % (self['run_name'], self['tag']))
                    local_dico['action'] += self.command_suggestion_html('pythia %s ' % self['run_name'])
            elif type == 'pythia':
                if self['tag'] == runresults.get_last_pythia():
                    if runresults.web:
                        local_dico['action'] = """
<FORM ACTION="http://%(web)s/cgi-bin/RunProcess/handle_runs-pl"  ENCTYPE="multipart/form-data" METHOD="POST">
<INPUT TYPE=HIDDEN NAME=directory VALUE="%(me_dir)s">
<INPUT TYPE=HIDDEN NAME=whattodo VALUE="remove_level">
<INPUT TYPE=HIDDEN NAME=level VALUE="pythia">
<INPUT TYPE=HIDDEN NAME=run VALUE="%(run_name)s">
<INPUT TYPE=HIDDEN NAME=tag VALUE=\"""" + self['tag'] + """\">
<INPUT TYPE=SUBMIT VALUE="Remove pythia">
</FORM>

<FORM ACTION="http://%(web)s/cgi-bin/RunProcess/handle_runs-pl"  ENCTYPE="multipart/form-data" METHOD="POST">
<INPUT TYPE=HIDDEN NAME=directory VALUE="%(me_dir)s">
<INPUT TYPE=HIDDEN NAME=whattodo VALUE="pgs">
<INPUT TYPE=HIDDEN NAME=run VALUE="%(run_name)s">
<INPUT TYPE=SUBMIT VALUE="Run Detector">
</FORM>"""
                    else:
                        local_dico['action'] = self.command_suggestion_html(
                                                'remove %s pythia --tag=%s' % \
                                                (self['run_name'], self['tag']))
                        local_dico['action'] += self.command_suggestion_html(
                         'pgs %(1)s or delphes %(1)s' % {'1': self['run_name']})
                else:
                    if runresults.web:
                        local_dico['action'] = ''
                    else:
                        local_dico['action'] = self.command_suggestion_html('remove %s  pythia --tag=%s'\
                                                                            % (self['run_name'], self['tag']))
            else:
                if runresults.web:
                    local_dico['action'] = """
<FORM ACTION="http://%(web)s/cgi-bin/RunProcess/handle_runs-pl"  ENCTYPE="multipart/form-data" METHOD="POST">
<INPUT TYPE=HIDDEN NAME=directory VALUE="%(me_dir)s">
<INPUT TYPE=HIDDEN NAME=whattodo VALUE="remove_level">
<INPUT TYPE=HIDDEN NAME=level VALUE=\"""" + str(type) + """\">
<INPUT TYPE=HIDDEN NAME=tag VALUE=\"""" + self['tag'] + """\">
<INPUT TYPE=HIDDEN NAME=run VALUE="%(run_name)s">
<INPUT TYPE=SUBMIT VALUE="Remove """ + str(type) + """\">
</FORM>"""
                else:
                    local_dico['action'] = self.command_suggestion_html('remove %s %s --tag=%s' %\
                                                              (self['run_name'], type, self['tag']))

            # create the text
            subresults_html += template % local_dico
            
            
        if subresults_html == '':
            if runresults.web:
                    action = """
<FORM ACTION="http://%(web)s/cgi-bin/RunProcess/handle_runs-pl"  ENCTYPE="multipart/form-data" METHOD="POST">
<INPUT TYPE=HIDDEN NAME=directory VALUE="%(me_dir)s">
<INPUT TYPE=HIDDEN NAME=whattodo VALUE="remove_level">
<INPUT TYPE=HIDDEN NAME=level VALUE="banner">
<INPUT TYPE=HIDDEN NAME=tag VALUE=\"""" + self['tag'] + """\">
<INPUT TYPE=HIDDEN NAME=run VALUE="%(run_name)s">
<INPUT TYPE=SUBMIT VALUE="Remove Banner">
</FORM>
                    
<FORM ACTION="http://%(web)s/cgi-bin/RunProcess/handle_runs-pl"  ENCTYPE="multipart/form-data" METHOD="POST">
<INPUT TYPE=HIDDEN NAME=directory VALUE="%(me_dir)s">
<INPUT TYPE=HIDDEN NAME=whattodo VALUE="banner">
<INPUT TYPE=HIDDEN NAME=run VALUE="%(run_name)s">
<INPUT TYPE=SUBMIT VALUE="Run the banner">
</FORM>"""
            else:
                    action = self.command_suggestion_html('remove %s banner --tag=%s' \
                                                                       % (self['run_name'], self['tag']))
                    action += self.command_suggestion_html('banner_run %s ' % self['run_name'])
            
            
            
            subresults_html = sub_part_template_parton % \
                          {'type': '', 
                           'run': self['run_name'],
                           'cross_span': 1,
                           'cross': self['cross'],
                           'err': self['error'],
                           'nb_event': self['nb_event'] and self['nb_event'] or 'No events yet',
                           'links': 'banner only',
                           'action': action
                           }                                
                                  
        if self.debug is KeyboardInterrupt:
            debug = '<br><font color=red>Interrupted</font>'
        elif isinstance(self.debug, basestring):
            if not os.path.isabs(self.debug) and not self.debug.startswith('./'):
                self.debug = './' + self.debug
            elif os.path.isabs(self.debug):
                self.debug = os.path.relpath(self.debug, self.me_dir)
            debug = '<br> <a href=\'%s\'> <font color=red>ERROR</font></a>' \
                                               % (self.debug)
        elif self.debug:
            text = str(self.debug).replace('. ','.<br>')
            if 'http' in text:
                pat = re.compile('(http[\S]*)')
                text = pat.sub(r'<a href=\1> here </a>', text)
            debug = '<br><font color=red>%s<BR>%s</font>' % \
                                           (self.debug.__class__.__name__, text)
        else:
            debug = ''                                       
        text = tag_template % {'tag_span': nb_line,
                           'run': self['run_name'], 'tag': self['tag'],
                           'subruns' : subresults_html,
                           'debug':debug}

        return text
        

    def command_suggestion_html(self, command):
        """return html button with code suggestion"""
        
        if command.startswith('pythia'):
            button = 'launch pythia'
        elif command.startswith('remove banner'):
            button = 'remove banner'
        elif command.startswith('remove'):
            button = 'remove run'
        elif command.startswith('banner_run'):
            button = 're-run from the banner'
        else:
            button = 'launch detector simulation'
        
        header = 'Launch ./bin/madevent in a shell, and run the following command: '
        return "<INPUT TYPE=SUBMIT VALUE='%s' onClick=\"alert('%s')\">" % (button, header + command)


        return  + '<br>'




