#!/usr/bin/env python


import gtk
import gtk.gdk as gdk
import pango
import gobject
import numpy as np
import random
import math
import os
import time
import pathfix
import atoms

PATH, POINT = range(2)
QUEUE_TYPE = 0
TYPE_ATOM, TYPE_LINE = range(2)
# For the following, DEPTH must be the next-to-last item.
ATOM_R, ATOM_RADIUS, ATOM_NUMBER, ATOM_ID, ATOM_DEPTH, ATOM_SIZE = range(1, 7)
LINE_R1, LINE_R2, LINE_COLOR, LINE_ID, LINE_WIDTH, LINE_DEPTH, LINE_SIZE = range(1, 8)
CONTEXT_ATOM, = range(1)


class atomview(gtk.Window):


#
# GUI -------------------------------------------------------------------------------------------
#


    def __init__(self):
        gtk.Window.__init__(self, gtk.WINDOW_TOPLEVEL)
        # Main window
        self.connect("destroy", self.event_close)
        self.connect("key_release_event", self.event_key_released)
        self.connect("key_press_event", self.event_key_pressed)
        self.set_resizable(True)
        self.set_default_size(402, 512)
        # Drawing area.
        self.area = gtk.DrawingArea()
        events = 0
        events = events | gdk.EXPOSURE_MASK | gdk.BUTTON_PRESS_MASK | gdk.BUTTON_RELEASE_MASK
        events = events | gdk.POINTER_MOTION_HINT_MASK | gdk.SCROLL
        self.area.set_events(events)                             
        self.area.connect("expose_event", self.event_exposed)
        self.area.connect("configure_event", self.event_configure)
        self.area.connect("button_press_event", self.event_button_press)
        self.area.connect("button_release_event", self.event_button_release)
        self.area.connect("motion_notify_event", self.event_mouse_move)
        self.area.connect("scroll_event", self.event_scroll)
        self.area.set_size_request(-1, -1)
        # VBox
        self.vbox = gtk.VBox()
        self.vbox.pack_start(self.area, True, True, 0)
        align = gtk.Alignment(0.5, 0.5)
        self.vbox.pack_start(align, False, False, 0)
        self.add(self.vbox)
        self.show_all()
        self.gui_members()
        self.event_configure()
        gobject.idle_add(self.gui_update)

    def gui_members(self):
        self.queue = []
        self.data = None
        self.lastr = None
        self.radius = 1.5
        self.scale = 8.0
        self.rotation = np.identity(3)
        self.button1 = False
        self.button2 = False
        self.button3 = False
        self.mouselast = (None, None)
        self.background = [1.0, 1.0, 1.0]
        self.keys = {}
        self.translate = np.array([0.0, 0.0, 16.0])
        self.gc = self.area.window.new_gc(foreground = gdk.Color(0, 0, 0))
        self.gcs = {}
        self.pixmap = None
        self.repeat = (1, 1, 1)
        self.boxsteps = 4
        self.maxfps = 60.0
        self.minspf = 1 / self.maxfps
        self.lastTime = time.time()
        self.highlights = {}
        self.bondscale = 3.0
        self.screenatoms = []
        self.grabbedatom = None
        self.mindepth = 2.0**16
        self.maxdepth = -2.0**16
        self.depthScale = 1.0
        
    def gui_update(self):
        wait = max(0, self.minspf - (time.time() - self.lastTime))      
        time.sleep(wait)                                                
        self.lastTime = time.time()       
        r = self.data.r
        if self.lastr == None:
            self.lastr = r
            self.event_exposed()
        elif r.shape == self.lastr.shape:
            if not (self.lastr == self.data.r).all():
                self.event_exposed()
        self.lastr = r
        return True

    
    def gui_key_on(self, key):
        return self.keys.has_key(key)


    def gui_mouse_atom(self):
        mx, my = self.mouselast
        mv = np.array([mx, my])
        atom = None
        index = 0
        for a in self.screenatoms:
            av = np.array([a[0], a[1]])
            dv = mv - av
            d = math.sqrt(np.dot(dv, dv))
            if d < a[2] and index < a[5]:
                atom = a[4]
                index = a[5]
        return atom
        

#
# EVENT -----------------------------------------------------------------------------------------
#


    def event_configure(self, *args):
        width, height = self.area.window.get_size()
        self.pixmap = gdk.Pixmap(self.area.window, width, height)      
        return True


    def event_key_pressed(self, widget, event):
        key = gdk.keyval_name(event.keyval)
        self.keys[key] = True
        if key == 'c' and (event.state & gdk.CONTROL_MASK):
            import sys
            sys.exit()


    def event_key_released(self, widget, event):
        key = gdk.keyval_name(event.keyval)
        if self.gui_key_on(key):
            self.keys.pop(key)


    def event_mouse_move(self, widget, event):
        mx, my, mask = self.area.window.get_pointer()
        if self.mouselast[0] == None:
            self.mouselast = (mx, my)
        else:
            dx = mx - self.mouselast[0]
            dy = my - self.mouselast[1]
            self.mouselast = (mx, my)
            if self.button1:
                if self.grabbedatom == None:
                    self.gfx_rot_x(dy * 0.0078125)
                    self.gfx_rot_y(dx * 0.0078125)
                else:
                    self.data[self.grabbedatom].r += np.dot(np.linalg.inv(self.rotation), np.array([dx, -dy, 0]) / self.scale / 2.0)
                self.event_exposed()
            elif self.button2:
                self.gfx_rot_z(-dx * 0.0078125)
                self.event_exposed()
            elif self.button3:
                self.translate += np.array([dx, -dy, 0]) / self.scale / 2
                self.event_exposed()
        print self.gui_mouse_atom()
        return True


    def event_button_press(self, widget, event):
        if event.button == 1:
            self.button1 = True
        if event.button == 2:
            self.button2 = True
        if event.button == 3:
            self.button3 = True
        if self.gui_key_on("m"):
            self.grabbedatom = self.gui_mouse_atom()
        return True


    def event_button_release(self, widget, event):
        if event.button == 1:
            self.button1 = False
        if event.button == 2:
            self.button2 = False
        if event.button == 3:
            self.button3 = False
        self.grabbedatom = None
        self.event_exposed()
        return True
            

    def event_scroll(self, widget, event):
        if self.gui_key_on("r"):
            if event.direction == gdk.SCROLL_UP:
                self.radius *= 1.1
                self.event_exposed()
            elif event.direction == gdk.SCROLL_DOWN:
                self.radius *= 0.9
                self.event_exposed()
        elif self.gui_key_on("f"):
            self.fadecheck.set_active(True)
            if event.direction == gdk.SCROLL_UP:
                self.depthScale *= 1.1
                self.event_exposed()
            elif event.direction == gdk.SCROLL_DOWN:
                self.depthScale *= 0.9
                self.event_exposed()
        elif self.gui_key_on("b"):
            if event.direction == gdk.SCROLL_UP:
                self.bondscale *= 1.1
            elif event.direction == gdk.SCROLL_DOWN:
                self.bondscale *= 0.9
            self.bondscale = min(self.bondscale, 0.75)
            self.event_exposed()
        else:
            if event.direction == gdk.SCROLL_UP:
                self.scale *= 1.1
            elif event.direction == gdk.SCROLL_DOWN:
                self.scale *= 0.9
            self.event_exposed()
        return True
                                                                                    

    def event_exposed(self, *args):
        self.gfx_clear()
        self.queue = []
        self.drawpoint = self.data
        self.gfx_queue_atoms()
        self.gfx_queue_box()
        self.gfx_transform_queue()
        self.gfx_sort_queue()
        self.gfx_draw_queue()
        self.gfx_draw_axes()
        self.area.window.draw_drawable(self.gc, self.pixmap, 0, 0, 0, 0, -1, -1)
        return True
                            

    def event_close(self, *args):
        gtk.main_quit()
        
        
#
# GRAPHICS --------------------------------------------------------------------------------------
#

    def gfx_set_color(self, r, g, b):
        self.gc.set_rgb_fg_color(gdk.Color(int(r * 65535), int(g * 65535), int(b * 65535)))
        
    def gfx_get_color_gc(self, r, g, b):
        rgb = (int(r * 65535), int(g * 65535), int(b * 65535))
        if rgb not in self.gcs:
            gc = self.area.window.new_gc()
            gc.set_rgb_fg_color(gdk.Color(rgb[0], rgb[1], rgb[2]))
            self.gcs[rgb] = gc
        return self.gcs[rgb]

    def gfx_set_line_width(self, width):
        self.gc.set_line_attributes(width, gdk.LINE_SOLID, gdk.CAP_ROUND, gdk.JOIN_ROUND)

    def gfx_queue_line(self, r1, r2, color, width = 1):
        line = [0] * LINE_SIZE
        line[QUEUE_TYPE] = TYPE_LINE
        line[LINE_R1] = np.copy(r1)
        line[LINE_R2] = np.copy(r2)
        line[LINE_COLOR] = color
        line[LINE_DEPTH] = (r1 + r2) / 2.0
        line[LINE_WIDTH] = width
        self.queue.append(line)
        
        
    def gfx_queue_atoms(self):
        r = self.drawpoint.r
        name = self.drawpoint.names
        # Try to get the box, if it exists.
        if self.data_have_box():
            box = self.drawpoint.box
        else:
            box = None
        # Queue the point
        for i in range(len(r)):
            atom = [0] * ATOM_SIZE
            atom[QUEUE_TYPE] = TYPE_ATOM
            atom[ATOM_R] = np.copy(r[i])
            atom[ATOM_RADIUS] = atoms.elements[name[i]]['radius']
            atom[ATOM_NUMBER] = atoms.elements[name[i]]['number']
            atom[ATOM_ID] = i % len(self.drawpoint)
            atom[ATOM_DEPTH] = 0
            self.queue.append(atom)
            
            
    def gfx_transform_queue(self):
        self.mindepth = 2.0**16
        self.maxdepth = -2.0**16
        r = self.drawpoint.r
        minx = min(r[:, 0])                         
        miny = min(r[:, 1])                         
        minz = min(r[:, 2])
        maxx = max(r[:, 0])
        maxy = max(r[:, 1])
        maxz = max(r[:, 2])
        midx = minx + (maxx - minx) / 2
        midy = miny + (maxy - miny) / 2
        midz = minz + (maxz - minz) / 2            
        mid = np.array([midx, midy, midz])
        r = self.drawpoint.r
        for i in range(len(self.queue)):
            q = self.queue[i]
            if q[QUEUE_TYPE] == TYPE_ATOM:
                q[ATOM_R] -= mid
                q[ATOM_R] = np.dot(self.rotation, q[ATOM_R])
                q[ATOM_R] += self.translate
                q[ATOM_DEPTH] = q[ATOM_R][2]
                self.mindepth = min(q[ATOM_DEPTH], self.mindepth)
                self.maxdepth = max(q[ATOM_DEPTH], self.maxdepth)
            else:                                   
                q[LINE_R1] -= mid
                q[LINE_R2] -= mid
                q[LINE_DEPTH] -= mid
                q[LINE_R1] = np.dot(self.rotation, q[LINE_R1])
                q[LINE_R2] = np.dot(self.rotation, q[LINE_R2])
                q[LINE_DEPTH] = np.dot(self.rotation, q[LINE_DEPTH])
                q[LINE_R1] += self.translate
                q[LINE_R2] += self.translate
                q[LINE_DEPTH] += self.translate
                q[LINE_DEPTH] = q[LINE_DEPTH][2]


    def gfx_sort_queue(self):
        def cmp_queue(a, b):
            if a[-1] > b[-1]:
                return 1
            else:
                return -1
        self.queue = sorted(self.queue, cmp_queue)
        

    def gfx_draw_queue(self):
        width, height = self.area.window.get_size()
        self.screenatoms = []
        index = 0
        for q in self.queue:
            index += 1
            if q[QUEUE_TYPE] == TYPE_ATOM:
                self.gfx_set_line_width(1)
                r = q[ATOM_R]
                rad = int(q[ATOM_RADIUS] * self.scale * self.radius)
                x = int(r[0] * self.scale * 2 + width * 0.5)
                y = int(-r[1] * self.scale * 2 + height * 0.5)
                if self.highlights.has_key(q[ATOM_ID]):
                    self.gfx_draw_circle(x, y, rad, q[ATOM_NUMBER], depth = q[ATOM_DEPTH], highlight = self.highlights[q[ATOM_ID]])
                else:
                    self.gfx_draw_circle(x, y, rad, q[ATOM_NUMBER], depth = q[ATOM_DEPTH], highlight = None)
                self.screenatoms.append([x, y, rad, q[ATOM_NUMBER], q[ATOM_ID], index])
            else:   
                q[LINE_R1][0] = q[LINE_R1][0] * self.scale * 2 + width * 0.5
                q[LINE_R1][1] = -q[LINE_R1][1] * self.scale * 2 + height * 0.5
                q[LINE_R2][0] = q[LINE_R2][0] * self.scale * 2 + width * 0.5
                q[LINE_R2][1] = -q[LINE_R2][1] * self.scale * 2 + height * 0.5
                self.gfx_draw_line(q[LINE_R1][0], q[LINE_R1][1], q[LINE_R2][0], q[LINE_R2][1], q[LINE_COLOR], q[LINE_WIDTH])


    def gfx_draw_axes(self):
        width, height = self.area.window.get_size()
        axes = np.identity(3) * 24
        axes[0] = np.dot(self.rotation, axes[0])
        axes[1] = np.dot(self.rotation, axes[1])
        axes[2] = np.dot(self.rotation, axes[2])
        x0 = 72
        y0 = height - 72
        self.gfx_draw_line(x0, y0, x0 - axes[0][0], y0 + axes[0][1], [0, 0, 0], 1)
        self.gfx_draw_line(x0, y0, x0 + axes[0][0], y0 - axes[0][1], [0, 0, 0], 1)
        self.gfx_draw_line(x0, y0, x0 - axes[1][0], y0 + axes[1][1], [0, 0, 0], 1)
        self.gfx_draw_line(x0, y0, x0 + axes[1][0], y0 - axes[1][1], [0, 0, 0], 1)
        self.gfx_draw_line(x0, y0, x0 - axes[2][0], y0 + axes[2][1], [0, 0, 0], 1)
        self.gfx_draw_line(x0, y0, x0 + axes[2][0], y0 - axes[2][1], [0, 0, 0], 1)
        font = pango.FontDescription("courier sans 10")
        X = self.area.create_pango_layout("+x")
        Y = self.area.create_pango_layout("+y")
        Z = self.area.create_pango_layout("+z")
        attr = pango.AttrList()
        fg = pango.AttrForeground(0, 0, 0, 0, -1)
        attr.insert(fg)
        X.set_alignment(pango.ALIGN_CENTER)
        Y.set_alignment(pango.ALIGN_CENTER)
        Z.set_alignment(pango.ALIGN_CENTER)
        X.set_font_description(font)
        Y.set_font_description(font)
        Z.set_font_description(font)
        X.set_attributes(attr)        
        Y.set_attributes(attr)        
        Z.set_attributes(attr)        
        self.gfx_set_color(self.background[0], self.background[1], self.background[2])
        self.pixmap.draw_layout(self.gc, int(x0 + axes[0][0]) - X.get_pixel_size()[0] / 2, 
                                int(y0 - axes[0][1]) - X.get_pixel_size()[1], X)
        self.pixmap.draw_layout(self.gc, int(x0 + axes[1][0]) - Y.get_pixel_size()[0] / 2, 
                                int(y0 - axes[1][1]) - Y.get_pixel_size()[1], Y)
        self.pixmap.draw_layout(self.gc, int(x0 + axes[2][0]) - Z.get_pixel_size()[0] / 2, 
                                int(y0 - axes[2][1]) - Z.get_pixel_size()[1], Z)

     
    def gfx_rot_x(self, theta):
        ct = math.cos(theta)
        st = math.sin(theta)
        m = np.array([[1, 0, 0], [0, ct, -st], [0, st, ct]])
        self.rotation = np.dot(m, self.rotation)
    
    
    def gfx_rot_y(self, theta):
        ct = math.cos(theta)
        st = math.sin(theta)
        m = np.array([[ct, 0, st], [0, 1, 0], [-st, 0, ct]])
        self.rotation = np.dot(m, self.rotation)
    
    
    def gfx_rot_z(self, theta):
        ct = math.cos(theta)
        st = math.sin(theta)
        m = np.array([[ct, -st, 0], [st, ct, 0], [0, 0, 1]])
        self.rotation = np.dot(m, self.rotation)
        

    def gfx_clear(self):
        self.gfx_set_color(self.background[0], self.background[1], self.background[2])
        width, height = self.area.window.get_size()
        self.pixmap.draw_rectangle(self.gc, True, 0, 0, width, height)


    def gfx_queue_box(self):
        try:
            self.drawpoint.box
        except:
            return
        bx = self.drawpoint.box
        b = np.array([[0, 0, 0],
                      [bx[1][0], bx[1][1], bx[1][2]],
                      [bx[1][0] + bx[0][0], bx[1][1] + bx[0][1], bx[1][2] + bx[0][2]],
                      [bx[0][0], bx[0][1], bx[0][2]],
                      [bx[2][0], bx[2][1], bx[2][2]],
                      [bx[2][0] + bx[1][0], bx[2][1] + bx[1][1], bx[2][2] + bx[1][2]],
                      [bx[2][0] + bx[1][0] + bx[0][0], bx[2][1] + bx[1][1] + bx[0][1], bx[2][2] + bx[1][2] + bx[0][2]],
                      [bx[2][0] + bx[0][0], bx[2][1] + bx[0][1], bx[2][2] + bx[0][2]]])
        index = [[0, 1], [0, 3], [0, 4], 
                 [7, 3], [7, 4], [7, 6], 
                 [5, 1], [5, 4], [5, 6],
                 [2, 6], [2, 3], [2, 1]]
        for i in index:
            r1 = b[i[0]]
            r2 = b[i[1]]
            for l in range(self.boxsteps):
                self.gfx_queue_line(r1 + (r2 - r1) * float(l) / self.boxsteps, 
                                    r1 + (r2 - r1) * float(l + 1) / self.boxsteps,
                                    [0, 0, 0])


    def gfx_draw_circle(self, x, y, r, element, depth = 1.0, highlight = None):
        depth = 1.0 - (depth - self.mindepth) / (self.maxdepth - self.mindepth)
        c = atoms.elements[element]['color']
        gc = self.gfx_get_color_gc(c[0] * 0.9, c[1] * 0.9, c[2] * 0.9)
        self.pixmap.draw_arc(gc, True, x - r, y - r, r * 2, r * 2, 0, 64 * 360)
        gc = self.gfx_get_color_gc(c[0], c[1], c[2])
        self.pixmap.draw_arc(gc, True, x - r / 2, y - r / 2, r / 2, r / 2, 0, 64 * 360)
        gc = self.gfx_get_color_gc(0, 0, 0)
        self.pixmap.draw_arc(gc, False, x - r, y - r, r * 2, r * 2, 0, 64 * 360)

    def gfx_draw_line(self, x1, y1, x2, y2, color, width = 1):
        self.gfx_set_color(color[0], color[1], color[2])
        self.gfx_set_line_width(width)
        self.pixmap.draw_line(self.gc, int(x1), int(y1), int(x2), int(y2))        


    def gfx_highlight(self, atomid, color):
        self.highlights[atomid] = color
        
    


#
# DATA ------------------------------------------------------------------------------------------
#


    def data_have_box(self):
        try:
            self.data.box
            return True
        except:
            return False

    def data_set(self, data):
        self.data = data
        self.event_exposed()





#
# MAIN ------------------------------------------------------------------------------------------
#


if __name__ == "__main__":
    import io
    import sys
    q = atomview()
    if len(sys.argv) > 1:
        q.data_set(io.loadcon(sys.argv[1]))
    gtk.main()



























