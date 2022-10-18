import numpy as np
from typing import List, Dict, Set,Tuple
import os
'''
oofempydir=''
#oofempydir='D:\\Development\\oofem\\build\\Debug'
oofempydir='D:\\Development\\oofem\\build\\Release'
if(oofempydir !=''):
    import sys
    sys.path.append(oofempydir)
try:
    import oofempy
except ImportError:
    pass
print('Used oofempy library: ' + oofempy.__file__)
'''
try:
    from femdir.oofemenum import *
except ImportError:
    from oofemenum import *

#



def get_referent_node_position(n1_p:np.ndarray, n2_p:np.ndarray, web_vec) ->np.ndarray:
    xvec = n2_p - n1_p
    y_vec = - np.cross(xvec,web_vec) # web is in the z direction
    y_vec = y_vec/np.linalg.norm(y_vec)
    ref_vec = (xvec+y_vec*np.linalg.norm(xvec))/2.0
    ref_node = n1_p + ref_vec
    return ref_node

def getDictIntVal(din: dict, key):
    val = din.get(key)
    if val is None:
        return val
    return int(val)

def getDictFloatVal(din: dict, key):
        val = din.get(key)
        if val is None:
            return val
        return float(val)
def get_range_from_items(items:List[int])->str:
    line = '{'
    sitems = sorted(items)
    first = sitems[0]
    last = sitems[0]
    for i in range(1,len(sitems)):
        cur= sitems[i]
        if cur > last+1:
            if last > first:
                line += ' ('+str(first)+' '+str(last)+')'
            else:
                line +=' '+str(first)
            first=cur
        last=cur
    if last > first:
        line += ' (' + str(first) + ' ' + str(last) + ')'
    else:
        line += ' ' + str(first)
    line += '}'
    return line
def get_list_from_items(items:List[int])->str:
    line = str(len(items))+' '
    for i in sorted(items):
        line += ' ' + str(i)
    return line

def get_list_adddata_from_items(items:List[int],adddata)->str:
    line = ''
    for i in sorted(items):
        line +=' '+str(i)+ ' ' + str(adddata)
    return line

class InputRecords():
    def __init__(self,numRecords=0):
        if numRecords == 0:
            self._lines=[]
        else:
            self._lines = [""]*numRecords

    @property
    def lines(self):
        return self._lines

    @staticmethod
    def getkeyname():
        return "warning: not implemented"

    def getNeutralFormatHelp(self):
        return None

    def getNeutralFormatRecordData(self,ixrec):
        return None

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        return False

    def getNumRecords(self):
        return len(self._lines)

    def getRecordLine(self,ixrec)->str:
        return  self._lines[ixrec]

    def setRecordLine(self,ixrec,line):
        self._lines[ixrec] = line

    def _tokenizeRecordLine(self, ixrec):
        return InputRecords.tokenize(self.getRecordLine(ixrec))

    def _dictinarizeRecordLine_11_1n(self, ixrec, n11keys)->(str, dict):
        sline = self._tokenizeRecordLine(ixrec)
        return sline[0].lower(), InputRecords._toDict1keyNvalue(sline, n11keys)

    def _dictinarizeRecordLine_11(self, ixrec) ->(str, dict):
        sline=self._tokenizeRecordLine(ixrec)
        return sline[0].lower(), InputRecords._toDict1key1value(sline)

    @staticmethod
    def tokenize(line):
        return line.split()

    @staticmethod
    def _toDict1key1value(sline):
        dout={}
        i=0
        while i < len(sline):
            dout[sline[i].lower()] = sline[i+1]
            i+=2
        return dout

    @staticmethod
    def _toDict1keyNvalue(sline, n11keys):
        idbasedkeywords = {'set','master','refnode','crosssect','ndofs','loadtype','cstype','outputatxy','outputtype',
                           'outputatz','outputcategory'}
        dout = {}
        i = 0
        while i < (n11keys * 2):
            dout[sline[i].lower()] = sline[i+1]
            i+=2
        while i < len(sline):
            if sline[i].lower() in idbasedkeywords:
                dout[sline[i].lower()] = sline[i + 1]
                i += 2
            else:
                n=int(sline[i+1])
                dout[sline[i].lower()] = sline[i+2:i+2+n]
                i += (n+2)
        return dout

class OutputManagerRecord(InputRecords):
    def __init__(self):
        super().__init__(1)

    @staticmethod
    def getkeyname():
        return 'out_manger'

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        #line='OutputManager tstep_all dofman_all element_all'
        line = 'OutputManager tstep_all dofman_all'
        #line = 'OutputManager'
        self.setRecordLine(0,line)
        return True




class AnalysisRecord(InputRecords):
    def __init__(self):
        super().__init__(1)
        #LinearStatic nsteps 3

    @staticmethod
    def getkeyname():
        return 'analysis'

    def getNeutralFormatRecordData(self, ixrec):
        sline = self._tokenizeRecordLine(ixrec)
        keyword = sline[0].lower()
        nlc = 0
        if  not(keyword == 'linearstatic'):
            print('Warning: record unhandled: ' + self.getRecordLine(ixrec))
            return keyword, 'Unhandled'
        for i in range(len(sline)):
            if sline[i] == 'nsteps':
                nlc = int(sline[i+1])
                break
        return keyword,nlc

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        (nlc) = data
        if nlc < 1:
            nlc = 1
        line='LinearStatic nsteps '+str(nlc)+' smtype 2 lstype 6'
        self.setRecordLine(0,line)
        return True

class DomainRecord(InputRecords):
    def __init__(self):
        super().__init__(1)

    @staticmethod
    def getkeyname():
        return 'domain'

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        line='domain 3dShell'
        self.setRecordLine(0,line)
        return True

class ComponentSizesRecord(InputRecords):
    def __init__(self):
        super().__init__(1)

    @staticmethod
    def getkeyname():
        return 'comp_sizes'


    def getNeutralFormatHelp(self):
        return ['num nodes', 'num elements', 'num props', 'num mat',
                'num bc and _loads', ' num bc cases', 'num groups','num init cond']

    def getNeutralFormatRecordData(self, ixrec):
        return self.getSizes()

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        line=''
        (ndofman, nelem, ncrosssect, nmat, nbc, nltf, nset, nic) = data
        line += 'ndofman '+ str(ndofman)+' nelem '+ str(nelem)+' ncrosssect '+ str(ncrosssect)
        line += ' nmat '+ str(nmat)+' nbc '+ str(nbc)+' nic '+ str(nic)
        line += ' nltf ' + str(nltf) + ' nset ' + str(nset)
        self.setRecordLine(0,line)
        return True

    def getSizes(self):
        keyword,dline = self._dictinarizeRecordLine_11(0)
        ndofman=nelem=ncrosssec=nmat=nbc=nic=nltf=nset=0
        val =getDictIntVal(dline,'ndofman')
        if val is not None:
            ndofman = val
        val = getDictIntVal(dline,'nelem')
        if val is not None:
            nelem = val
        val = getDictIntVal(dline,'ncrosssect')
        if val is not None:
            ncrosssec = val
        val = getDictIntVal(dline,'nmat')
        if val is not None:
            nmat = val
        val = getDictIntVal(dline,'nbc')
        if val is not None:
            nbc = val
        val = getDictIntVal(dline,'nic')
        if val is not None:
            nic = val
        val = getDictIntVal(dline,'nltf')
        if val is not None:
            nltf = val
        val = getDictIntVal(dline,'nset')
        if val is not None:
            nset = val
        return ndofman, nelem, ncrosssec, nmat, nbc, nltf, nset,nic

class NodeRecords(InputRecords):
    def __init__(self, numRecords):
        super().__init__(numRecords)

    @staticmethod
    def getkeyname():
        return 'nodes'

    def getNeutralFormatRecordData(self, ixrec):
        keyword, dline = self._dictinarizeRecordLine_11_1n(ixrec, 1)
        id = -1
        coords = []
        lcs=[]
        haveLCS = False
        if keyword == 'node':
            for key, value in dline.items():
                if key == keyword:
                    id = int(value)
                elif key == 'coords':
                    for x in value:
                        coords.append(float(x))
                elif key == 'lcs':
                    haveLCS = True
                    for x in value:
                        lcs.append(float(x))
            if haveLCS:
                return keyword, id,coords,lcs
            else:
                return keyword, id, coords
        elif keyword == 'rigidarmnode':
            for key, value in dline.items():
                if key == keyword:
                    id = int(value)
                elif key == 'coords':
                    for x in value:
                        coords.append(float(x))
                elif key == 'master':
                    idmaster = int(value)
                elif key == 'dofidmask':
                    dofidmask = []
                    for item in value:
                        dofidmask.append(int(item))
                elif key == 'doftype':
                    doftype = []
                    for item in value:
                        doftype.append(int(item))
                elif key == 'mastermask':
                    mastermask = []
                    for item in value:
                        mastermask.append(int(item))
            return keyword, id, coords,idmaster,dofidmask # doftype and mastermask are determined by dofidmask
        else:
            print('Warning: record unhandled: ' + self.getRecordLine(ixrec))
            return keyword,'Unhandled'

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        if data[0] == 'node':

            keyword= data[0]
            id = data[1]
            coords = data[2]
            line=keyword
            line+=' '+str(id)
            line += ' coords '+str(len(coords))
            for item in coords:
                line +=' '+str(item)
            if len(data) > 3:
                dofidmask = data[3]
                line += ' ndofs ' + str(len(dofidmask))
                line += ' dofidmask ' + str(len(dofidmask))
                for item in dofidmask:
                    line += ' ' + str(item)
            self.setRecordLine(ixrec,line)
            return True
        elif data[0] == 'rigidarmnode':
            (keyword, id, coords,idmaster,dofidmask) = data
            line = keyword
            line += ' ' + str(id)
            line += ' coords ' + str(len(coords))
            for item in coords:
                line += ' ' + str(item)
            line += ' master ' + str(idmaster)
            line += ' dofidmask ' + str(len(dofidmask))
            for item in dofidmask:
                line += ' ' + str(item)
            line += ' doftype ' + str(len(dofidmask))
            for item in dofidmask:
                line += ' 2'
            line += ' mastermask ' + str(len(dofidmask))
            for item in dofidmask:
                line += ' 1'
            self.setRecordLine(ixrec, line)
            return True
        else:
            print('Warning! Node keyword unhandled:' + str(data))

class ElementRecords(InputRecords):
    def __init__(self, numRecords):
        super().__init__(numRecords)

    @staticmethod
    def getkeyname():
        return 'elements'
    
    def getNeutralFormatRecordData(self, ixrec):
        keyword, dline = self._dictinarizeRecordLine_11_1n(ixrec, 1)
        idprop=-1
        id=-1
        nodes = []
        for key, value in dline.items():
            if key == keyword:
                id = int(value)
            elif key == 'nodes':
                nn=len(value)
                for i in range(nn):
                    nodes.append(int(value[i]))
            elif key == 'refnode':
                nodes.append(int(value))
            elif key == 'crosssect':
                idprop = int(value)
        if idprop > -1:
            return keyword, id,nodes,idprop
        else:
            return keyword, id, nodes
    def setRecordDataNeutralFormat(self,ixrec,data:()):
        if len(data)== 4:
            (keyword, id,nodes,idprop) = data
        else:
            (keyword, id, nodes) = data
        line = keyword
        line += ' '+ str(id)
        if str(keyword).lower().startswith('beam'):
            line += ' nodes ' + str(len(nodes)-1)
            for item in nodes[:-1]:
                line += ' ' + str(item)
            line += ' refnode ' + str(nodes[-1])
        else:
            line += ' nodes ' + str(len(nodes))
            for item in nodes:
                line += ' ' + str(item)
        if len(data) == 4:
            line += ' crosssect ' + str(idprop)
        line += get_nip_str_for_element(keyword)
        self.setRecordLine(ixrec, line)
        return True

class SetRecords(InputRecords):
    def __init__(self, numRecords):
        super().__init__(numRecords)

    @staticmethod
    def getkeyname():
        return 'sets'

    def getNeutralFormatRecordData(self, ixrec):
        sline = self._tokenizeRecordLine(ixrec)
        id=int(sline[1])
        idlist = []
        type = sline[2].lower()
        if type == 'nodes' or type == 'elements':
            for idit in range(4,len(sline)):
                idlist.append(int(sline[idit]))
        elif type == 'noderanges' or type == 'elementranges':
            #{(1 6) 10 11 15 16 20 21 25}
            line = self.getRecordLine(ixrec)
            icbs = line.find('{')
            icbe = line.find('}')
            sline= InputRecords.tokenize(line[icbs+1:icbe])
            i=0
            while i < len(sline):
                if sline[i].find('(') == 0:
                    istart= int(sline[i][1:])
                    iend = int(sline[i+1][:-1])
                    for idit in range(istart,iend+1):
                        idlist.append(idit)
                    i+=2
                else:
                    idlist.append(int(sline[i]))
                    i+=1
        elif type == 'elementboundaries':
            type = 'elements'
            npairs= int(int(sline[3])/2)
            for idit in range(npairs):
                idlist.append(int(sline[(idit*2)+4]))
        elif type == 'allelements':
            pass
        else:
            print('Warning! set type unhandled:' + type)
            return type, 'Unknown'

        if type == 'noderanges':
            type = 'nodes'
        elif type == 'elementranges':
            type = 'elements'
        return type,id,idlist

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        (type, id, idlist) = data
        maxitems = 10
        line='set ' +str(id)+' '
        if len(idlist) > maxitems and type != 'pressureelements':
            type=type[:-1]+'ranges'
        if type == 'pressureelements':
            type = 'elementboundaries'
            line+= type+' '+str(len(idlist)*2)
            line += get_list_adddata_from_items(idlist,1)
        elif type == 'allelements':
            line += type
        else:
            line += type + ' '

            if len(idlist) > maxitems:
                line+= get_range_from_items(idlist)
            else:
                line+= get_list_from_items(idlist)
        self.setRecordLine(ixrec,line)
        return True

class CrossSectionRecords(InputRecords):
    def __init__(self, numRecords):
        super().__init__(numRecords)

    @staticmethod
    def getkeyname():
        return 'cross_sects'

    def getNeutralFormatRecordData(self, ixrec):
        cstype, dline = self._dictinarizeRecordLine_11(ixrec)
        null = None
        proptype = null

        if cstype.lower() == 'simplecs':
            id = getDictIntVal(dline,cstype)
            idmat = getDictIntVal(dline, 'material')
            idset = getDictIntVal(dline, 'set')
            tp = getDictFloatVal(dline, 'thick')
            if tp is not null:

                add_char_candidates = {'reldrillstiffness','drilltype'}
                add_chars={}
                for key in add_char_candidates:
                    value = getDictFloatVal(dline, key)
                    if value is not null:
                        add_chars[key] = value
                if(len(add_chars)>0):
                    proptype = 'plate_add'  # beam prop with stiffness parameters only
                    return proptype, id, idmat, tp, idset,add_chars
                else:
                    proptype = 'plate'  # beam prop with stiffness parameters only
                    return proptype, id, idmat, tp, idset
            area = getDictFloatVal(dline, 'area')
            if area is not null:
                Iy = getDictFloatVal(dline, 'iy')
                Iz = getDictFloatVal(dline, 'iz')
                Ik = getDictFloatVal(dline, 'ik')
                if (Iy is not null) or (Iz is not null) or (Ik is not null):
                    proptype = 'k_beam' # beam prop with stiffness parameters only
                    shearareay = getDictFloatVal(dline, 'shearareay')
                    shearareaz = getDictFloatVal(dline, 'shearareaz')
                    shear_coeff = getDictFloatVal(dline, 'beamshearcoeff')
                    sfiff_list = [area,Iy,Iz,Ik,shear_coeff,shearareay,shearareaz]
                    for i in range(len(sfiff_list)):
                        if sfiff_list[i] is null:
                            sfiff_list[i]=0.0
                    return proptype,id,idmat,sfiff_list,idset
                else:
                    proptype = 'k_rod'  # rod prop with stiffness parameters only
                    return proptype, id, idmat, area, idset
        print('Warning: record unhandled: ' + self.getRecordLine(ixrec))
        return None, 'Unhandled'

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        line = 'simplecs '
        proptype=data[0]
        if proptype =='plate':
            (proptype, id, idmat, tp, idset) = data
            line += str(id)
            line += ' thick '+str(tp)
            line += ' material ' + str(idmat)
            if idset is not None:
                line += ' set ' + str(idset)
        elif proptype =='plate_add':
            (proptype, id, idmat, tp, idset,dict_add_cschar) = data
            line += str(id)
            line += ' thick '+str(tp)
            line += ' material ' + str(idmat)
            if idset is not None:
                line += ' set ' + str(idset)
            for key,value in dict_add_cschar.items():
                    line += ' '+key+' ' + str(value)
        elif proptype =='k_beam':
            (proptype,id,idmat,sfiff_list,idset) = data
            line += str(id)
            for key,stiff in sfiff_list.items():
                if stiff > 0.0:
                    line += ' '+key+' ' + str(stiff)
            line += ' material ' + str(idmat)
            if idset is not None:
                line += ' set ' + str(idset)
        elif proptype =='k_rod':
            (proptype,id,idmat,area,idset) = data
            line += str(id)
            line += ' area ' + str(area)
            line += ' material ' + str(idmat)
            if idset is not None:
                line += ' set ' + str(idset)
        else:
            print('Warning: proptype '+ str(proptype) + ' unhandled: ' + str(data))
            return False
        self.setRecordLine(ixrec,line)
        return True

class MaterialRecords(InputRecords):
    def __init__(self, numRecords):
        super().__init__(numRecords)

    @staticmethod
    def getkeyname():
        return 'materials'

    def getNeutralFormatRecordData(self, ixrec):
        mattype, dline = self._dictinarizeRecordLine_11(ixrec)
        if mattype.lower() == 'isole':
            id = getDictIntVal(dline,mattype)
            E = getDictFloatVal(dline,'e')
            ni = getDictFloatVal(dline, 'n')
            rho = getDictFloatVal(dline, 'd')
            return mattype,id,E,ni,rho
        return None

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        (mattype, id, E, ni, rho) = data
        line=mattype +' '+str(id)+' E '+str(E)+' n '+str(ni)+' d '+str(rho)+' talpha 0'
        self.setRecordLine(ixrec,line)
        return True

class BCandLoadRecords(InputRecords):
    def __init__(self, numRecords):
        super().__init__(numRecords)

    @staticmethod
    def getkeyname():
        return 'bcs_loads'

    def getNeutralFormatRecordData(self, ixrec):
        keyword, dline = self._dictinarizeRecordLine_11_1n(ixrec, 2)
        id = int(dline[keyword])
        idfun = dline.get('loadtimefunction')
        if keyword == 'nodalload' or keyword == 'boundarycondition':
            lvals = []
            ldofs = []
            idset = int(dline.get('set'))
            if idfun is not None:
                idfun = int(idfun)
                sindices = dline['dofs']
                if keyword == 'nodalload':
                    scomponents = dline['components']
                elif keyword == 'boundarycondition':
                    scomponents = dline['values']
                else:
                    print('Warning: record unhandled: ' + self.getRecordLine(ixrec))
                for i in range(len(sindices)):
                    lvals.append(float(scomponents[i]))
                    ldofs.append(int(sindices[i]))
                return keyword, id, idset, idfun, ldofs, lvals
        elif keyword == 'deadweight':
            lvals = []
            idset = int(dline.get('set'))
            if idfun is not None:
                idfun = int(idfun)
                scomponents = dline['components']
                for i in range(len(scomponents)):
                    lvals.append(float(scomponents[i]))
                return keyword, id, idset, idfun, lvals
        elif keyword == 'constantsurfaceload':
            lvals = []
            idset = int(dline.get('set'))
            if idfun is not None:
                idfun = int(idfun)
                loadtype = int(dline['loadtype'])
                cstype = int(dline['cstype'])
                scomponents = dline['components']
                is_pressure = True
                press = float(scomponents[2])  # only z (z axis is in a direction of plate element normal)
                for i in range(len(scomponents)):
                    # only pressure load is handled
                    if i != 2 and float(scomponents[i]) != 0.0:
                        print('Warning: components unhandled: ' + scomponents)
                        is_pressure = False
                        break
                if cstype != 1:
                    print('Warning: cstype unhandled: ' + dline['cstype'])
                    is_pressure = False
                if loadtype != 3:
                    print('Warning: loadtype unhandled: ' + dline['loadtype'])
                    is_pressure = False
                if is_pressure:
                    return keyword, id, idset, idfun, press
        print ('Warning: record unhandled: ' + self.getRecordLine(ixrec))
        return keyword, 'Unhandled'

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        keyword = data[0].lower()
        id = data[1]
        idset = data[2]
        idfun = data[3]
        if keyword == 'pressureload':
            keyword = 'constantsurfaceload'
        if keyword == 'accelerationload':
            keyword = 'deadweight'
        line = keyword+ ' ' + str(id) + ' loadTimeFunction ' + str(idfun)

        if keyword ==  'nodalload':
            dofs = data[4]
            vals = data[5]
            line += ' dofs ' + str(len(dofs))
            for dof in dofs:
                line += ' ' + str(dof)
            line+= ' components ' + str(len(dofs))
            for val in vals:
                line += ' ' + str(val)
            line += ' csType 0'
        elif keyword ==  'boundarycondition':
            dofs = data[4]
            vals = data[5]
            line += ' dofs ' + str(len(dofs))
            for dof in dofs:
                line += ' ' + str(dof)
            line+= ' values ' + str(len(dofs))
            for val in vals:
                line += ' ' + str(val)
        elif keyword == 'deadweight':
            vals = data[4]
            line += ' components ' + str(len(vals))
            for val in vals:
                line += ' ' + str(val)
        elif keyword ==  'constantsurfaceload':
            press = data[4]
            line += ' loadType 3 ndofs 6 components 6 0.0 0.0 {} 0.0 0.0 0.0 csType 1'.format(press)
        else:
            print('Warning! load and Boundary Condition keyword unhandled:' + str(data))
            return False
        line += ' set ' + str(idset)
        self.setRecordLine(ixrec, line)
        return True


class InitialConditionRecords(InputRecords):
    def __init__(self, numRecords):
        super().__init__(numRecords)

    @staticmethod
    def getkeyname():
        return 'init_conds'

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        line=''
        self.setRecordLine(0,line)
        return True

class TimeFunctionRecords(InputRecords):
    def __init__(self, numRecords):
        super().__init__(numRecords)

    @staticmethod
    def getkeyname():
        return 'tfuncs'

    def getNeutralFormatRecordData(self, ixrec):
        sline = self._tokenizeRecordLine(ixrec)
        id = int(sline[1])
        if sline[0].lower()=='constantfunction':
            fval = float(sline[3])
            idlc = -999
            if np.round(fval,1) != 1.0:
                print('Warning f(t) not 1.0: ' + self.getRecordLine(ixrec))
        elif sline[0].lower()=='piecewiselinfunction':
            keyword, dline = self._dictinarizeRecordLine_11_1n(ixrec, 2)
            n = int(dline['npoints'])
            idlc = -1
            for i in range(n):
                fval = np.round(float(dline['f(t)'][i]),1)
                if np.round(fval,1) == 1.0:
                    idlc = int(dline['t'][i])
        return id,idlc

    def setRecordDataNeutralFormat(self,ixrec,data:()):
        (id,idlc,nlc) = data
        if idlc == -999:
            line = 'ConstantFunction '+ str(id) +' f(t) 1.0'
            self.setRecordLine(ixrec, line)
        else:
            line = 'piecewiselinfunction ' + str(id) + ' nPoints '+str(nlc)+' t '+str(nlc)
            for i in range(nlc):
                line+= ' '+str(i+1)
            line+= ' f(t) '+str(nlc)
            for i in range(nlc):
                if i == (idlc-1):
                    line+= ' 1'
                else:
                    line+= ' 0'


        self.setRecordLine(ixrec, line)
        return True

class OOFEMInputFile():
    def __init__(self,fileName):
        if fileName is not None:
            self.filename = fileName
        else:
            self.filename = ""
        self._outfilename= ""
        self._description =""
        self._records = {}
        if self.filename != "":
            self._readtorecords()

    def getRecordGroup(self,key)->InputRecords:
        return self._records.get(key)

    def _addRecordGroup(self, record_group)->InputRecords:
        self._records[record_group.getkeyname()] = record_group
        return record_group



    def _readtorecords(self):
        self._records.clear()
        self._outfilename = ""
        self._description = ""
        with open(self.filename, 'r') as f:
            lines = f.readlines()
            shift =0
            while (lines[shift][0]=="#" or lines[shift].strip() ==""):
                shift+=1
            self._outfilename = lines[shift].strip()
            shift += 1
            while (lines[shift][0]=="#" or lines[shift].strip() ==""):
                shift+=1
            self._description = lines[shift]
            shift += 1
            # analysis
            rg = self._addRecordGroup(AnalysisRecord())
            shift = self._addLinesToRecordGroup(rg, 1, shift, lines)
            # errorcheck
            while (lines[shift][0]=="#" or lines[shift].strip() ==""):
                shift+=1
            if lines[shift].strip().startswith('errorcheck'):
                shift += 1
            # domain
            rg = self._addRecordGroup(DomainRecord())
            shift = self._addLinesToRecordGroup(rg, 1, shift, lines)
            # out_manager
            rg = self._addRecordGroup(OutputManagerRecord())
            shift = self._addLinesToRecordGroup(rg, 1, shift, lines)
            # sizes
            rg = self._addRecordGroup(ComponentSizesRecord())
            shift = self._addLinesToRecordGroup(rg, 1, shift, lines)
            ndofman, nelem, ncrosssec, nmat, nbc, nltf, nset,nic = rg.getSizes()
            # nodes
            rg = self._addRecordGroup(NodeRecords(ndofman))
            shift = self._addLinesToRecordGroup(rg, ndofman, shift, lines)
            #elements
            rg = self._addRecordGroup(ElementRecords(nelem))
            shift = self._addLinesToRecordGroup(rg, nelem, shift, lines)
            # support for set after elements
            while (lines[shift][0]=="#" or lines[shift].strip() ==""):
                shift+=1
            if nset > 0 and lines[shift].strip().lower().startswith('set'):
                rg = self._addRecordGroup(SetRecords(nset))
                shift = self._addLinesToRecordGroup(rg, nset, shift, lines)
                nset = 0
            # cross section
            rg = self._addRecordGroup(CrossSectionRecords(ncrosssec))
            shift = self._addLinesToRecordGroup(rg, ncrosssec, shift, lines)
            # materials
            rg = self._addRecordGroup(MaterialRecords(nmat))
            shift = self._addLinesToRecordGroup(rg, nmat, shift, lines)
            # bc & _loads
            rg = self._addRecordGroup(BCandLoadRecords(nbc))
            shift = self._addLinesToRecordGroup(rg, nbc, shift, lines)
            # initial conditions
            rg = self._addRecordGroup(InitialConditionRecords(nic))
            shift = self._addLinesToRecordGroup(rg, nic, shift, lines)
            # time functions
            rg = self._addRecordGroup(TimeFunctionRecords(nltf))
            shift = self._addLinesToRecordGroup(rg, nltf, shift, lines)
            if nset > 0:
                rg = self._addRecordGroup(SetRecords(nset))
                self._addLinesToRecordGroup(rg, nset, shift, lines)

    def init_input_records(self,data:()):
        """
        Initialize all input records grops using size data (the same as necessary for ComponentSizesRecord)
        The function also  prepares DomainRecord, OutputManagerRecord and ComponentSizesRecord

        :param data: (ndofman, nelem, ncrosssect, nmat, nbc, nltf, nset, nic) = data
        :return:
        """
        (ndofman, nelem, ncrosssect, nmat, nbc, nltf, nset, nic) = data
        rg = self._addRecordGroup(AnalysisRecord())
        # domain
        rg = self._addRecordGroup(DomainRecord())
        rg.setRecordDataNeutralFormat(0, (0))  # no input necessary
        # out_manager
        rg = self._addRecordGroup(OutputManagerRecord())
        rg.setRecordDataNeutralFormat(0, (0)) # no input necessary
        # sizes
        rg = self._addRecordGroup(ComponentSizesRecord())
        rg.setRecordDataNeutralFormat(0,data)
        # nodes
        rg = self._addRecordGroup(NodeRecords(ndofman))
        # elements
        rg = self._addRecordGroup(ElementRecords(nelem))
        # cross section
        rg = self._addRecordGroup(CrossSectionRecords(ncrosssect))
        # materials
        rg = self._addRecordGroup(MaterialRecords(nmat))
        # bc & _loads
        rg = self._addRecordGroup(BCandLoadRecords(nbc))
        # initial conditions
        rg = self._addRecordGroup(InitialConditionRecords(nic))
        # time functions
        rg = self._addRecordGroup(TimeFunctionRecords(nltf))
        #sets
        rg = self._addRecordGroup(SetRecords(nset))

    def _check_input_records(self):
        error = 0
        msg = ''
        null = None
        # identify pressure groups
        press_groups = set()
        drg = self.getRecordGroup(BCandLoadRecords.getkeyname())  # data record group
        if drg is not null:
            n = drg.getNumRecords()
            for i in range(n):
                neutral = drg.getNeutralFormatRecordData(i)
                idfun = neutral[3]
                keyword = neutral[0]
                if keyword == 'constantsurfaceload':
                    keyword, id, idset, idfun, press = neutral
                    press_groups.add(idset)
        if len(press_groups) > 0:
            # Sets
            drg = self.getRecordGroup(SetRecords.getkeyname())  # data record group
            n = drg.getNumRecords()
            for i in range(n):
                sline = drg._tokenizeRecordLine(i)
                id = int(sline[1])
                type = sline[2].lower()
                if id in press_groups and type != 'elementboundaries':
                    type,id,idlist = drg.getNeutralFormatRecordData(i)
                    type = 'elementboundaries'
                    line = 'set '+ str(id)+ ' '  +type + ' ' + str(len(idlist) * 2)
                    line += get_list_adddata_from_items(idlist, 1)
                    drg.setRecordLine(i, line)
                    msg+='\nSet {} transformed to elementboundaries'.format(id)


        return error,msg
    def write_from_records(self,file_path):
        error, msg = self._check_input_records()
        filenameext = os.path.basename(file_path)
        filename, file_extension = os.path.splitext(filenameext)
        outfilename = filename+'.out'
        description = "Automaticaly generated file for  static linear analysis of ship structure"
        lines=[]
        lines.append(outfilename+'\n')
        lines.append(description+'\n')
        # analysis
        self._addRecordGroupToLines(self.getRecordGroup(AnalysisRecord.getkeyname()),lines)
        # domain
        self._addRecordGroupToLines(self.getRecordGroup(DomainRecord.getkeyname()), lines)
        # out_manager
        self._addRecordGroupToLines(self.getRecordGroup(OutputManagerRecord.getkeyname()), lines)
        # sizes
        self._addRecordGroupToLines(self.getRecordGroup(ComponentSizesRecord.getkeyname()), lines)
        # nodes
        self._addRecordGroupToLines(self.getRecordGroup(NodeRecords.getkeyname()), lines)
        #elements
        self._addRecordGroupToLines(self.getRecordGroup(ElementRecords.getkeyname()), lines)
        # cross section
        self._addRecordGroupToLines(self.getRecordGroup(CrossSectionRecords.getkeyname()), lines)
        # materials
        self._addRecordGroupToLines(self.getRecordGroup(MaterialRecords.getkeyname()), lines)
        # bc & _loads
        self._addRecordGroupToLines(self.getRecordGroup(BCandLoadRecords.getkeyname()), lines)
        # initial conditions
        self._addRecordGroupToLines(self.getRecordGroup(InitialConditionRecords.getkeyname()), lines)
        # time functions
        self._addRecordGroupToLines(self.getRecordGroup(TimeFunctionRecords.getkeyname()), lines)
        # sets
        self._addRecordGroupToLines(self.getRecordGroup(SetRecords.getkeyname()), lines)
        with open(file_path, 'w') as f:
            f.writelines(lines)
        return error,msg


    def _addLinesToRecordGroup(self, rg, nlines, shift, lines):
        i=0
        ir=0
        while i < nlines:
            if (lines[i + shift][0]=='#' or lines[i+shift].strip() == ""):
                nlines+=1
            else:
                rg.setRecordLine(ir, lines[i + shift])
                ir += 1
            i+=1
        return shift+nlines

    def _addRecordGroupToLines(self, rg:InputRecords, lines:List[str]):
        if rg is not None:
            for line in rg.lines:
                lines.append(line+'\n')
'''
class OOFEMAnalysisModel():
    def __init__(self, fileName,dict_idset_outtypes:Dict[int,OutputElementType]=None):
        null = None
        if fileName is null:
            self._file_path = ''
        self._eng_model:oofempy.EngngModel = null
        self._vtkxml_nodes:oofempy.VTKXMLExportModule = null
        self._vtkxml_shell:oofempy.VTKXMLExportModule = null
        self._vtkxml_beam:oofempy.VTKXMLExportModule = null
        self._nlc:int = 0
        self._file_path,self._dr,self._eng_model,self._nlc = self.init_problem(fileName)
        self._wdir = os.getcwd()
        self._dict_idset_outtypes:Dict[int,OutputElementType] = dict_idset_outtypes


    @property
    def num_lc(self):
        return self._nlc

    def get_cs(self,id_cs:int):
        domain = self._eng_model.giveDomain(1)
        cs = domain.giveCrossSection(id_cs)
        return cs

    def init_problem(self,file_path)->(oofempy.EngngModel,oofempy.VTKXMLExportModule,oofempy.VTKXMLExportModule,oofempy.VTKXMLExportModule):
        null = None
        try:
            if not os.path.exists(file_path):
                print('Error, file path does not exist: ' + file_path)
                return file_path,null, null,0
            dr = oofempy.OOFEMTXTDataReader(file_path)
            eng_mdl: oofempy.EngngModel = oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, False, None, False)
            nlc = eng_mdl.giveNumberOfSteps()
        except BaseException as err:
            print(f"Unexpected {err=}, {type(err)=}")
            file_path=None
            dr=None
            eng_mdl = None
            nlc =0
        return file_path, dr, eng_mdl,nlc

    def _check_num_elements_in_set(self,domain,id_set)->int:
        set = domain.giveSet(id_set)
        lst = set.giveElementList()
        #have_elements = lst.isEmpty()
        num_elements = lst.giveSize()
        return num_elements


    def init_analyisis(self, dict_idset_outtypes:Dict[int,OutputElementType]):
        list_shell=[]
        list_beam = []
        domain = self._eng_model.giveDomain(1)
        nset = domain.giveNumberOfSets()  # sets can not have unused ids!
        if dict_idset_outtypes is not None:
            for idset, outtype in dict_idset_outtypes.items():
                if idset > nset:
                    continue
                if idset < 0: # additional option to handle the case where output sets are added at the end of existing sets
                    idset = nset + 1 - idset
                if self._check_num_elements_in_set(domain,idset) == 0:
                    continue
                if outtype is OutputElementType.Shell:
                    list_shell.append(idset)
                elif outtype is OutputElementType.Beam:
                    list_beam.append(idset)

        self._wdir = os.getcwd()
        os.chdir(os.path.dirname(self._file_path))

        self._vtkxml_nodes = oofempy.vtkxml(1, self._eng_model, domain_all=True, tstep_all=True, dofman_all=True,
                                            element_all=True, primvars=(1,), pythonExport=1)

        if len(list_shell) == 0:
            self._vtkxml_shell = None
        else:
            tuple_shell = tuple(list_shell)
            self._vtkxml_shell = oofempy.vtkxml(2, self._eng_model, domain_all=True, tstep_all=True, dofman_all=True,
                                                element_all=False, cellvars=(200,), stype=1,
                                                regionsets=tuple_shell, pythonExport=1)

        if len(list_beam) == 0:
            self._vtkxml_beam = None
        else:
            tuple_beam = tuple(list_beam)
            self._vtkxml_beam = oofempy.vtkxml(3, self._eng_model, domain_all=True, tstep_all=True, dofman_all=True,
                                               element_all=False, cellvars=(201,),
                                               regionsets=tuple_beam, pythonExport=1)

    def finalize_analysis(self):
        os.chdir(self._wdir)

    def analyse_next_loadcase(self):
        null = None
        self._eng_model.preInitializeNextStep()
        self._eng_model.giveNextStep()
        currLC = self._eng_model.giveCurrentStep()
        self._eng_model.initializeYourself(currLC)
        self._eng_model.solveYourselfAt(currLC)
        self._eng_model.updateYourself(currLC)
        self._eng_model.terminate(currLC)
        disp_node_out = self._vtkxml_nodes.getPrimaryVars()['DisplacementVector']
        if self._vtkxml_shell is null:
            shell_el_out = null
        else:
            shell_el_out = self._vtkxml_shell.getCellVars()['IST_Shell_SxSyTxy_Top_Bottom']
        if self._vtkxml_beam is null:
            beam_el_out = null
        else:
            beam_el_out = self._vtkxml_beam.getCellVars()['IST_Beam_FxMyMz_Start_End']
        reactions_dict = self._vtkxml_nodes.getReactionForces()
        reactions = [reactions_dict['nodeid'],
                     reactions_dict['dofid'],
                     reactions_dict['reaction']]
        nodes_def_data=self._vtkxml_nodes.getNodes()['1']
        node_ids = []
        for nod_data in nodes_def_data:
            node_ids.append(nod_data[2])
        return disp_node_out, shell_el_out, beam_el_out,reactions,node_ids


    def analyse_all_loadcases(self, dict_idset_outtypes:Dict[int,OutputElementType]):
        node_out=[]
        shell_out = []
        beam_out = []
        react_out = []
        node_id_out = []
        try:
            self.init_analyisis(dict_idset_outtypes)
            for ilc in range(self.num_lc):
                disp_node_out, shell_el_out, beam_el_out, reactions,node_ids =  self.analyse_next_loadcase()
                node_out.append(disp_node_out)
                shell_out.append(shell_el_out)
                beam_out.append(beam_el_out)
                react_out.append(reactions)
                node_id_out.append(node_ids)
            self.finalize_analysis()
        except BaseException as err:
            print(f"Unexpected {err=}, {type(err)=}")

        return node_out, shell_out, beam_out, react_out, node_id_out



    def analyse(self):
        return self.analyse_all_loadcases(self._dict_idset_outtypes)
'''