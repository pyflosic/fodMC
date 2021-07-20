from ase.io import read, write 
import numpy as np
from mayavi import mlab

from traits.api import Any, HasTraits, Button, Instance, Range, Str 
from traitsui.api import View, Group, Item,TextEditor

#import sys 
#import io
#from contextlib import redirect_stdout
import os 
try: 
    import fodmc 
except: print('You have not given the path to fodmc.so') 

# author: Sebastian Schwalbe 
# date:   07.01.2020 
# description: gui for atomic systems to define bonds 

# Notes: (1) Currenly the bonds are shifted along z problem if molecule is oriented in z direction 
#        (2) Lone electrons currenlty "randomly" arrange around a atom 
#            this needs to be updated that it corresponds to Lewis like lone electrons/pairs 
#        (3) pygui.redo function need be recoded, may collect all bonds and lone pairs as objects 
#            and then delete the objects by identifiers. Currently the redo function 
#            deletes both the last bond and the last lone electron 
#        (4) we need to set up more elements in pygui (see settings) 

# How to use gui elements with mayavi 
# /home/schwalbe/.local/lib/python3.6/site-packages/mayavi/tools/animator.py 
# traits texteditor 
# https://github.com/enthought/traitsui/blob/master/examples/demo/Standard_Editors/TextEditor_demo.py

class pygui_output(HasTraits):
    # This a gui element handling out but from pygui 
    output = Str("      ") 
    
    traits_view = View(Group(Item('output',style='custom',show_label=False),Item('_'),label=''),title='Info')
    
    def __init(self):
        HasTraits.__init__(self)
        self.ui = None

    def show(self):
        self.ui = self.configure_traits()

    def close(self):
        if self.ui is not None:
          self.ui.dispose()

class pygui_utils(HasTraits):
    # This is a gui element providing buttons/events 
    info =  Button('Info')
    get = Button('Get Bonding')
    run = Button('Run fodmc') 
    redo = Button('Redo')
    
    traits_view = View(Group(
                       Item('info'),
                       Item('get'),
                       Item('run'),
                       Item('redo'),
                       show_labels=False),
                       Item('_'),
                       title='PyGUI',
                       buttons=['OK'])

    def __init__(self,pygui):
        HasTraits.__init__(self)
        self.pygui = pygui 
        self.ui = None 

    def show(self):
        self.ui = self.edit_traits()

    def close(self):
      if self.ui is not None:
          self.ui.dispose()

    def _get_fired(self):
        # write and read fodmc input 
        # from pygui bonds/lone_elec 
        self.write_fodmc()
        # system -- fodmc input file 
        f = open('system','r') 
        output = f.read() 
        f.close()
        out = pygui_output(output=output)
        out.show()

    def _redo_fired(self):
        output = str(self.pygui.redo())

    def _info_fired(self):
        # Some information for pygui 
        s1 = 'Selection: action\n'
        s2 = '1 atom: click -> adds one electron\n'
        s3 = '2 atoms: 1st -> single bond\n'
        s4 = '2 atoms: 2nd -> double bond\n'
        s5 = '2 atoms: 3rd -> trible bond' 
        output = s1 + s2 + s3 + s4 + s5  
        out = pygui_output(output=output)
        out.show()

    def _run_fired(self):
        try:
            # magic to capture that output:
            # from http://stackoverflow.com/questions/977840/redirecting-fortran-called-via-f2py-output-in-python
            #      http://websrv.cs.umt.edu/isis/index.php/F2py_example
            output_file = 'fodmc.out'
            # open outputfile
            outfile = os.open(output_file, os.O_RDWR|os.O_CREAT)
            # save the current file descriptor
            save = os.dup(1)
            # put outfile on 1
            os.dup2(outfile, 1)
            # end magic 
            
            # Fortran call
            fodmc.fodmc_mod.get_guess()
#            fodmc.fodmc.points_on_sphere_metropolis_spin_centers()
            #
            
            # restore the standard output file descriptor
            os.dup2(save, 1)
            # close the output file
            os.close(outfile)
            f = open(output_file,'r') 
            output = f.read()
        except: output = 'You have not provided the path to fodmc.so!'
        out = pygui_output(output=output)
        out.show()

    def write_fodmc(self):
        # the actual fodmc input write function 
        bond_settings = {'single' : '(1-1)', 'double' : '(2-2)' , 'trible': '(3-3)'}
        lone_elec_settings =  {'single' : '(1-0)', 'double' : '(1-1)' , 'trible': '(2-1)'}
        # system : fodmc input 
        o = open('system','w') 
        natoms = len(self.pygui.atoms)
        sys = self.pygui.atoms.get_chemical_formula()
        sym = self.pygui.atoms.get_chemical_symbols()
        pos = self.pygui.atoms.get_positions() 
        o.write('%i %s \n' %(natoms,sys))
        # old 
        o.write('angstrom\n')
        # new 
        #o.write('angstrom fix1s\n')
        for s in range(len(sym)):
            o.write('%s %0.9f %0.9f %0.9f\n' %(sym[s],pos[s][0],pos[s][1],pos[s][2]))
        o.write('con_mat\n') 
        bonds = self.pygui.bonds
        b_str = ''
        for b in list(bonds.keys()): 
            idx1 = bonds[b]['start_index']+1
            idx2 = bonds[b]['end_index']+1
            t = bonds[b]['type']
            b_str = b_str + '(%i-%i)-%s ' %(idx1,idx2,bond_settings[t])
        o.write(b_str+'\n')
        lone_elec = self.pygui.lone_elec 
        l_str = ''
        for l in list(lone_elec):
            idx1 = lone_elec[l]['start_index']+1
            t = lone_elec[l]['type']
            l_str = l_str + '%i-%s ' %(idx1,lone_elec_settings[t])
        o.write(l_str+'\n')
        o.write('\n')
        o.close() 

class pygui:
    # a grafical inferface for bonding visualization 

    # scale_factor : radius of the sphere visualized 
    settings = {'O' : {'scale_factor' : 1.2, 'color': (1,0,0)},
                'H' : {'scale_factor' : 0.3, 'color': (1,1,1)},
                'C' : {'scale_factor' : 1.1, 'color': (0.5,0.5,0.5)},
                'N' : {'scale_factor' : 1.15, 'color': (0,0.75,0.75)},
                'He' : {'scale_factor' : 0.3, 'color': (0,0.5,0)},
                'X' : {'scale_factor' : 0.3, 'color': (0.5,0,0)}
                }   
                

    def __init__(self,f_struct):
        self.f_struct = f_struct 
        # read struct 
        self.atoms = read(f_struct)
        
        # mayavi obj 
        self.initialize() 

    def initialize(self):

        self.figure = mlab.figure('PyGUI', bgcolor=(0, 0, 0), size=(350, 350),)
        f = mlab.gcf()
        self.f = f 
        #mlab.clf()
        self.figure.scene.disable_render = True

        atoms_glyph = []
        a = 0
        for a in self.atoms: 
            # Debug 
            # print(a.position)
            # print(a.symbol)
            tmp_atom = self.make_atom(a) 
            atoms_glyph.append(tmp_atom) 
        self.atoms_glyph = atoms_glyph 

        # Add an outline to show the selected point and center it on the first
        # data point.
        # outline = mlab.outline(line_width=3)
        # outline.outline_mode = 'cornered'
        # outline.bounds = (self.atoms[0].position[0]-0.1, self.atoms[0].position[0]+0.1,
        #                   self.atoms[0].position[1]-0.1, self.atoms[0].position[1]+0.1,
        #                   self.atoms[0].position[2]-0.1, self.atoms[0].position[2]+0.1)
        # self.outline = outline 
        
        # Every object has been created, we can reenable the rendering.
        self.figure.scene.disable_render = False

        # Here, we grab the points describing the individual glyph, to figure
        # out how many points are in an individual glyph.
        self.glyph_points = len(self.atoms_glyph) #.glyph.glyph_source.glyph_source.output.points.to_array()

        self.selected_atoms = []
        self.selected_atoms_glyph = []
        self.bonds = {}
        # running index for bonds 
        self.bi = 0
        self.lone_elec = {}
        # running index for lone electrons
        self.li = 0
        self.picker = self.figure.on_mouse_pick(self.picker_callback)
        # A keypress method which might useful later 
        #self.f.scene.interactor.add_observer('KeyPressEvent', self.key_callback)

        # Decrease the tolerance, so that we can more easily select a precise
        # point.
        self.picker.tolerance = 0.01

        # Title visualized in pygui 
        mlab.title('Click')

        # addition gui elements using pygui_utils class 
        ui = pygui_utils(self)
        self.ui = ui 

    def make_atom(self,atom):
        # draw a atom as a sphere using mlab 
        atom = mlab.points3d(atom.position[0], atom.position[1], atom.position[2],
                  scale_factor=self.settings[atom.symbol]['scale_factor'],
                  resolution=20,
                  color=self.settings[atom.symbol]['color'],
                  scale_mode='none')
        return atom 

    def redo(self):
        # redo n-1 action 
        # this is work in progress and needs to be optimized 

        lone_elec = self.lone_elec
        lone_elec_keys = list(self.lone_elec.keys())
        
        # Debug 
        # print('lone_elec_keys',lone_elec_keys)
        # print(self.li) 
        try:
            lone_elec_keys.remove('e_%s'%(self.li-1))
            lone_elec.pop('e_%s'%(self.li-1),'None')
            li = int(lone_elec_keys[-1].split('_')[-1])
        except: 'Nothing' #li = self.li+1  
        # Debug 
        # print('lone_elec_keys',lone_elec_keys)
        picker = self.picker 
        #lone_elec = self.lone_elec
        bonds = self.bonds 
        bond_keys = list(self.bonds.keys())
        # Debug 
        # print('bond_keys',bond_keys)
        # print(self.bi)
        try:
            bond_keys.remove('b_%s'%(self.bi-1))
            bonds.pop('b_%s'%(self.bi-1),None)
            bi = int(bond_keys[-1].split('_')[-1])
        except: 'Nothing' #bi = self.bi+1 
        # Debug 
        # print('----------------------------------------------------------')
        # print('b_%s'%(self.bi-1))
        # print(bond_keys)
        mlab.close()
        bi = self.bi 
        li = self.li 
        self.initialize()
        self.bi = bi 
        self.li = li
        #self.__init__(self.f_struct)
        self.picker = picker 
        self.lone_elec = lone_elec 
        self.bonds = bonds 
        if len(lone_elec_keys) >= 1:
            for l in list(lone_elec_keys):
                pos1 = self.lone_elec[l]['start']
                pos2 = self.lone_elec[l]['end']
                offset = self.lone_elec[l]['offset']
                color = self.lone_elec[l]['color']
                x = np.array([pos1[0],pos2[0]])
                y = np.array([pos1[1],pos2[1]])
                z = np.array([pos1[2],pos2[2]])
                # single lone elec 
                if self.lone_elec[l]['type'] == 'single':
                    mlab.points3d(x-offset, y-offset,z+offset,
                                  scale_factor=0.1,
                                  resolution=20,
                                  color=color,
                                  scale_mode='none')
                if self.lone_elec[l]['type'] == 'double':
                    mlab.points3d(x-offset, y-offset,z+offset,
                                  scale_factor=0.1,
                                  resolution=20,
                                  color=color,
                                  scale_mode='none')
                    mlab.points3d(x+offset, y+offset,z+offset,
                                  scale_factor=0.1,
                                  resolution=20,
                                  color=color,
                                  scale_mode='none')
                if self.lone_elec[l]['type'] == 'triple':
                    mlab.points3d(x-offset, y-offset,z+offset,
                                  scale_factor=0.1,
                                  resolution=20,
                                  color=color,
                                  scale_mode='none')
                    mlab.points3d(x+offset, y+offset,z+offset,
                                  scale_factor=0.1,
                                  resolution=20,
                                  color=color,
                                  scale_mode='none')
                    mlab.points3d(x+offset, y+offset,z-offset,
                                  scale_factor=0.1,
                                  resolution=20,
                                  color=color,
                                  scale_mode='none')
                if self.lone_elec[l]['type'] == 'quad':
                    mlab.points3d(x-offset, y-offset,z+offset,
                                  scale_factor=0.1,
                                  resolution=20,
                                  color=color,
                                  scale_mode='none')
                    mlab.points3d(x+offset, y+offset,z+offset,
                                  scale_factor=0.1,
                                  resolution=20,
                                  color=color,
                                  scale_mode='none')
                    mlab.points3d(x+offset, y+offset,z-offset,
                                  scale_factor=0.1,
                                  resolution=20,
                                  color=color,
                                  scale_mode='none')
                    mlab.points3d(x-offset, y-offset,z+offset,
                                  scale_factor=0.1,
                                  resolution=20,
                                  color=color,
                                  scale_mode='none')

        if len(bond_keys) >=1:
            for b in list(bond_keys):
                pos1 = self.bonds[b]['start']
                pos2 = self.bonds[b]['end']
                x = np.array([pos1[0],pos2[0]])
                y = np.array([pos1[1],pos2[1]])
                z = np.array([pos1[2],pos2[2]])

                # single bound
                if self.bonds[b]['type'] == 'single':
                    mlab.plot3d(x,y,z,tube_radius=0.025, colormap='Spectral')
                # double bound
                if self.bonds[b]['type'] == 'double':
                    mlab.plot3d(x,y,z,tube_radius=0.025, colormap='Spectral')
                    mlab.plot3d(x,y,z+0.1,tube_radius=0.025, colormap='Spectral')
                # trible bound 
                if self.bonds[b]['type'] == 'trible':
                    mlab.plot3d(x,y,z,tube_radius=0.025, colormap='Spectral')
                    mlab.plot3d(x,y,z+0.1,tube_radius=0.025, colormap='Spectral')
                    mlab.plot3d(x,y,z-0.1,tube_radius=0.025, colormap='Spectral')

    def key_callback(self,vtk_obj,event):
         # a key_press method might useful later 
         print(vtk_obj.GetKeyCode())

    def picker_callback(self,picker):
        # selecting one atom  -> adding lone electron 
        # selecting two atoms -> adding bonds 
       
        # Hide the title 
        # Title visualized in pygui
        mlab.title('')

        for a in range(len(self.atoms_glyph)):
            if picker.actor in self.atoms_glyph[a].actor.actors:
                # yes we picked a point 
                # Debug 
                # print('yes')
                x, y, z = self.atoms[a].position[0], self.atoms[a].position[1], self.atoms[a].position[2]
                # the selected point is colored green (0,1,0) 
                self.atoms_glyph[a].actor.property.color = (0,1,0)
                # position which is selected 
                # Debug 
                # print('pos',x,y,z)
                # Move the outline to the data point.
                #self.outline.bounds = (x-0.1, x+0.1,y-0.1, y+0.1,z-0.1, z+0.1)
                # We collect the selected poibts as ase.atoms and 
                # as mlab object 
                self.selected_atoms.append(self.atoms[a])
                self.selected_atoms_glyph.append(self.atoms_glyph[a])

        # if we have selected two times 
        if len(self.selected_atoms) == 2:
            pos1 = self.selected_atoms[0].position
            pos2 = self.selected_atoms[1].position
            x = np.array([pos1[0],pos2[0]])
            y = np.array([pos1[1],pos2[1]])
            z = np.array([pos1[2],pos2[2]])
            
            # lone electrons 
            # if we have selected two times the same atom 
            if (pos1 == pos2).all():
                offset = 0.3*self.settings[self.selected_atoms[0].symbol]['scale_factor']
                mlab.points3d(self.selected_atoms[0].position[0]-offset, self.selected_atoms[0].position[1]-offset,self.selected_atoms[0].position[2]+offset,
                              scale_factor=0.1,
                              resolution=20,
                              color=self.settings[self.selected_atoms[0].symbol]['color'],
                              scale_mode='none')
                t = 'single'
                update = True 
                del_idx = []
                for l in list(self.lone_elec.keys()):
                    if (self.lone_elec[l]['start'] == pos1).all() or (self.lone_elec[l]['start'] == pos2).all():
                        if (self.lone_elec[l]['end'] == pos2).all() or (self.lone_elec[l]['end'] == pos1).all():
                            # double elec
                            if self.lone_elec[l]['type'] == 'single':
                                del_idx.append(l)
                                mlab.points3d(self.selected_atoms[0].position[0]+offset, 
                                        self.selected_atoms[0].position[1]+offset,
                                        self.selected_atoms[0].position[2]+offset,
                                        scale_factor=0.1,
                                        resolution=20,
                                        color=self.settings[self.selected_atoms[0].symbol]['color'],
                                        scale_mode='none')
                                t = 'double'
                                update = True 
                            # trible elec
                            if self.lone_elec[l]['type'] == 'double':
                                del_idx.append(l)
                                mlab.points3d(self.selected_atoms[0].position[0]+offset, 
                                        self.selected_atoms[0].position[1]+offset,
                                        self.selected_atoms[0].position[2]-offset,
                                        scale_factor=0.1,
                                        resolution=20,
                                        color=self.settings[self.selected_atoms[0].symbol]['color'],
                                        scale_mode='none')
                                t = 'trible'
                                update = True 
                            # quad elec
                            if self.lone_elec[l]['type'] == 'trible':
                                del_idx.append(l)
                                mlab.points3d(self.selected_atoms[0].position[0]-offset,
                                        self.selected_atoms[0].position[1]-offset,
                                        self.selected_atoms[0].position[2]-offset,
                                        scale_factor=0.1,
                                        resolution=20,
                                        color=self.settings[self.selected_atoms[0].symbol]['color'],
                                        scale_mode='none')
                                t = 'quad'
                                update = True 
                            #if self.lone_elec[l]['type'] == 'quad':
                            #    self.redo() 
                            #    #self.picker = self.figure.on_mouse_pick(self.picker_callback)
                            #    #self.picker.tolerance = 0.01
                            #    update = False 
                            #    #self.show()
                            #    #return
                                
                if update == True:
                    # delete all update bonds
                    for l in del_idx:
                        del self.lone_elec[l]
                    lone_elec = {'e_%s' %(self.li):{'type':t,
                                                    'start':self.selected_atoms[0].position,
                                                    'end':self.selected_atoms[1].position,
                                                    'start_index':self.selected_atoms[0].index,
                                                    'end_index':self.selected_atoms[1].index,
                                                    'offset':offset,
                                                    'color':self.settings[self.selected_atoms[0].symbol]['color']
                                                    }}
                    self.li = self.li +1
                    self.lone_elec.update(lone_elec)
                    # Debug 
                    # print(self.lone_elec)
            # bonds 
            # if we have selected two different atoms 
            if (pos1 != pos2).any():
                # Debug 
                # print('bond')
                # single bond 
                mlab.plot3d(x,y,z,tube_radius=0.025, colormap='Spectral')
                t = 'single' 
                update = True
                del_idx = []
                for b in list(self.bonds.keys()):
                    if (self.bonds[b]['start'] == pos1).all() or (self.bonds[b]['start'] == pos2).all():
                        if (self.bonds[b]['end'] == pos2).all() or (self.bonds[b]['end'] == pos1).all():
                            # double bound 
                            if self.bonds[b]['type'] == 'single':
                                del_idx.append(b) 
                                mlab.plot3d(x,y,z+0.1,tube_radius=0.025, colormap='Spectral')
                                t = 'double' 
                                update = True 
                            # trible bound 
                            if self.bonds[b]['type'] == 'double':
                                del_idx.append(b)
                                mlab.plot3d(x,y,z-0.1,tube_radius=0.025, colormap='Spectral')
                                t = 'trible'
                                update = True 
                            #if self.bonds[b]['type'] == 'trible':
                            #    self.redo()
                            #    #self.picker = self.figure.on_mouse_pick(self.picker_callback)
                            #    #self.picker.tolerance = 0.01
                            #    update = False 
                            #    print(dir(self))
                            #    #self.picker_callback
                            #    #self.show()
                            #    #return 
                
                if update == True:
                    # delete all update bonds 
                    for k in del_idx:
                        del self.bonds[k]

                    bond = {'b_%s' %(self.bi):{'type':t,
                                               'start':self.selected_atoms[0].position,
                                               'end':self.selected_atoms[1].position,
                                               'start_index' :self.selected_atoms[0].index,
                                               'end_index' :self.selected_atoms[1].index
                                               }}
                    self.bi = self.bi +1 
                    self.bonds.update(bond)
                    # Debug 
                    # print(self.bonds)
    
        if len(self.selected_atoms) >= 3:
                # reset colors if we have more then two selected atoms 
                for a in range(len(self.selected_atoms)):
                    self.selected_atoms_glyph[a].actor.property.color = self.settings[self.selected_atoms[a].symbol]['color']
                self.selected_atoms = []
                self.selected_atoms_glyph = []
        self.ui.close()    
        self.ui.show()

    def show(self):
        # main gui 
        mlab.show() 

if __name__ == "__main__":
    #f_struct = 'H2O.xyz'
    #f_struct = 'CH4.xyz' 
    #f_struct = 'C6H6.xyz'
    # fodmc bonding test 
    #f_struct = 'C2H2.xyz'
    # fodmc lone elec test 
    #f_struct = 'NH3.xyz'
    f_struct = 'Nuc_FOD.xyz'
    pg =pygui(f_struct=f_struct)
    pg.show()


