#!/usr/bin/env python

import os
import math
import time

import gtk
import gtk.gdk as gdk
import gtk.glade as glade
import gobject

import numpy as np

import pathfix
import atoms


class queueitem:
    def __init__(self, kind):
        self.kind = kind

class atomview(gtk.Window):
#
# GUI -------------------------------------------------------------------------------------------
#

    def __init__(self, gui = None):
        gtk.Window.__init__(self, gtk.WINDOW_TOPLEVEL)
        # Main window
        self.connect("destroy", self.event_close)
        self.connect("key_release_event", self.event_key_released)
        self.connect("key_press_event", self.event_key_pressed)
        self.set_resizable(True)
        # Glade
        gladetree = gtk.glade.XML(os.path.join(pathfix.path, "tools/atomview.glade"))
        gladewindow = gladetree.get_widget("window")
        self.moviescale = gladetree.get_widget("moviescale")
        self.playbutton = gladetree.get_widget("playbutton")
        self.pausebutton = gladetree.get_widget("pausebutton")
        self.boxbutton = gladetree.get_widget("boxbutton")
        self.resetbutton = gladetree.get_widget("resetbutton")
        self.fps = gladetree.get_widget("fps")
        self.holder = gladetree.get_widget("holder")
        if gui is not None:
            self.holder.pack_start(gui, False, False, 0)
        # Events
        self.playbutton.connect("clicked", self.event_toggle_play, True)
        self.pausebutton.connect("clicked", self.event_toggle_play, False)
        self.resetbutton.connect("clicked", lambda w: self.gfx_reset_transform())
        self.boxbutton.connect("clicked", lambda w: self.queue_draw())
        # Drawing area.
        self.area = gladetree.get_widget("atomview")
        self.area.connect("expose_event", self.event_exposed)
        self.area.connect("configure_event", self.event_configure)
        self.area.connect("button_press_event", self.event_button_press)
        self.area.connect("button_release_event", self.event_button_release)
        self.area.connect("motion_notify_event", self.event_mouse_move)
        self.area.connect("scroll_event", self.event_scroll)
        self.area.set_size_request(512, 512)
        self.add(gladewindow)
        self.show_all()
        self.gui_members()
        self.event_configure()
        self.gfx_setup_colors()
        self.gfx_reset_transform()
        # Timeout
        gobject.timeout_add(8, self.event_timeout)

    def gui_members(self):
        self.last_draw = 0.0
        self.playing = False
        self.queue = []
        self.data = None
        self.button1 = False
        self.button2 = False
        self.button3 = False
        self.mouselast = (None, None)
        self.background = [0.949, 0.945, 0.941]
        self.keys = {}
        self.black_gc = self.area.get_style().black_gc
        self.white_gc = self.area.get_style().white_gc
        self.background_gc = self.gfx_get_color_gc(self.background[0], self.background[1], self.background[2])
        self.pixmap = None
        self.repeat = (1, 1, 1)
        self.lastTime = time.time()
        self.screenatoms = []
        self.colors = []

    def gui_key_on(self, key):
        return key in self.keys

#
# EVENT -----------------------------------------------------------------------------------------
#


    def event_configure(self, *args):
        x, y, self.width, self.height = self.area.get_allocation()
        self.pixmap = gdk.Pixmap(self.area.window, self.width, self.height)
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
                self.gfx_rot_x(dy * 0.009)
                self.gfx_rot_y(dx * 0.009)
                self.queue_draw()
            elif self.button2:
                self.gfx_rot_z(-dx * 0.0078125)
                self.queue_draw()
            elif self.button3:
                self.translate += np.array([dx, -dy, 0]) / self.scale / 2
                self.queue_draw()
        return True


    def event_button_press(self, widget, event):
        if event.button == 1:
            self.button1 = True
        if event.button == 2:
            self.button2 = True
        if event.button == 3:
            self.button3 = True
        return True


    def event_button_release(self, widget, event):
        if event.button == 1:
            self.button1 = False
        if event.button == 2:
            self.button2 = False
        if event.button == 3:
            self.button3 = False
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
        self.drawpoint = self.data[0]
        if len(self.data) > 1:
            self.drawpoint = self.data[int(self.moviescale.get_value())]
        self.gfx_queue_atoms()
        if self.boxbutton.get_active():
            self.gfx_queue_box()
        self.gfx_transform_queue()
        self.gfx_sort_queue()
        self.gfx_draw_queue()
        self.area.window.draw_drawable(self.white_gc, self.pixmap, 0, 0, 0, 0, self.width, self.height)
        self.last_draw = time.time()

    def event_toggle_play(self, widget, play):
        self.playing = play

    def event_timeout(self):
        if self.playing and len(self.data) > 1:
            if time.time() - self.last_draw > 1.0/self.fps.get_value():
                nextframe = self.moviescale.get_value() + 1
                if nextframe > len(self.data) - 1:
                    nextframe = 0
                self.moviescale.set_value(nextframe)
                self.queue_draw()
        return True

    def event_close(self, *args):
        gtk.main_quit()


#
# GRAPHICS --------------------------------------------------------------------------------------
#

    def gfx_get_color_gc(self, r, g, b):
        rgb = (int(r * 65535), int(g * 65535), int(b * 65535))
        gc = self.area.window.new_gc()
        gc.set_rgb_fg_color(gdk.Color(rgb[0], rgb[1], rgb[2]))
        return gc

    def gfx_setup_colors(self):
        for i in range(atoms.numElements):
            c = atoms.elements[i]['color']
            self.colors.append(self.gfx_get_color_gc(c[0], c[1], c[2]))

    def gfx_reset_transform(self):
        self.radius = 1.5
        self.scale = 8.0
        self.rotation = np.identity(3)
        self.translate = np.array([0.0, 0.0, 16.0])
        self.queue_draw()

    def gfx_queue_line(self, r1, r2, color, width = 1):
        line = queueitem("line")
        line.r1 = np.copy(r1)
        line.r2 = np.copy(r2)
        line.color = color
        line.depth = (r1 + r2) / 2.0
        line.width = width
        self.queue.append(line)

    def gfx_queue_atoms(self):
        r = self.drawpoint.r
        name = self.drawpoint.names
        for i in range(len(r)):
            atom = queueitem("atom")
            atom.r = np.copy(r[i])
            atom.radius = atoms.elements[name[i]]['radius']
            atom.number = atoms.elements[name[i]]['number']
            atom.id = i % len(self.drawpoint)
            atom.depth = 0
            self.queue.append(atom)

    def gfx_queue_box(self):
        try:
            self.drawpoint.box
        except:
            return
        bx = self.drawpoint.box
        b = np.array([[0, 0, 0], [bx[1][0], bx[1][1], bx[1][2]], [bx[1][0] +
                     bx[0][0], bx[1][1] + bx[0][1], bx[1][2] + bx[0][2]],
                     [bx[0][0], bx[0][1], bx[0][2]], [bx[2][0], bx[2][1],
                     bx[2][2]], [bx[2][0] + bx[1][0], bx[2][1] + bx[1][1],
                     bx[2][2] + bx[1][2]], [bx[2][0] + bx[1][0] + bx[0][0],
                     bx[2][1] + bx[1][1] + bx[0][1], bx[2][2] + bx[1][2] +
                     bx[0][2]], [bx[2][0] + bx[0][0], bx[2][1] + bx[0][1],
                     bx[2][2] + bx[0][2]]])
        index = [[0, 1], [0, 3], [0, 4], [7, 3], [7, 4], [7, 6], [5, 1], [5, 4],
                 [5, 6], [2, 6], [2, 3], [2, 1]]
        for i in index:
            r1 = b[i[0]]
            r2 = b[i[1]]
            boxsteps = 4
            for l in range(boxsteps):
                self.gfx_queue_line(r1 + (r2 - r1) * float(l) / boxsteps,
                                    r1 + (r2 - r1) * float(l + 1) /
                                    boxsteps, [0, 0, 0])

    def gfx_center_atoms(self):
        r = self.data[0].r
        minx = min(r[:, 0])
        miny = min(r[:, 1])
        minz = min(r[:, 2])
        maxx = max(r[:, 0])
        maxy = max(r[:, 1])
        maxz = max(r[:, 2])
        midx = minx + (maxx - minx) / 2
        midy = miny + (maxy - miny) / 2
        midz = minz + (maxz - minz) / 2
        self.center = np.array([midx, midy, midz])


    def gfx_transform_queue(self):
        r = self.drawpoint.r
        for i in range(len(self.queue)):
            q = self.queue[i]
            if q.kind == "atom":
                q.r -= self.center
                q.r = np.dot(self.rotation, q.r)
                q.r += self.translate
                q.depth = q.r[2]
            else:
                q.r1 -= self.center
                q.r2 -= self.center
                q.depth -= self.center
                q.r1 = np.dot(self.rotation, q.r1)
                q.r2 = np.dot(self.rotation, q.r2)
                q.depth = np.dot(self.rotation, q.depth)
                q.r1 += self.translate
                q.r2 += self.translate
                q.depth += self.translate
                q.depth = q.depth[2]


    def gfx_sort_queue(self):
        def cmp_queue(a, b):
            if a.depth > b.depth:
                return 1
            else:
                return -1
        self.queue = sorted(self.queue, cmp_queue)


    def gfx_draw_queue(self):
        s2 = self.scale * 2
        w2 = self.width * 0.5
        h2 = self.height * 0.5
        for q in self.queue:
            if q.kind == "atom":
                r = q.r
                rad = int(q.radius * self.scale * self.radius)
                x = int(r[0] * s2 + w2)
                y = int(-r[1] * s2 + h2)
                self.gfx_draw_circle(x, y, rad, q.number)
            else:
                q.r1[0] = q.r1[0] * self.scale * 2 + self.width * 0.5
                q.r1[1] = -q.r1[1] * self.scale * 2 + self.height * 0.5
                q.r2[0] = q.r2[0] * self.scale * 2 + self.width * 0.5
                q.r2[1] = -q.r2[1] * self.scale * 2 + self.height * 0.5
                self.gfx_draw_line(q.r1[0], q.r1[1], q.r2[0], q.r2[1])


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
        self.pixmap.draw_rectangle(self.background_gc, True, 0, 0, self.width, self.height)

    def gfx_draw_circle(self, x, y, r, element):
        r = max(1,r)
        self.pixmap.draw_arc(self.colors[element], True, x - r, y - r, r * 2, r * 2, 0, 64 * 360)
        self.pixmap.draw_arc(self.black_gc, False, x - r, y - r, r * 2, r * 2, 0, 64 * 360)

    def gfx_draw_line(self, x1, y1, x2, y2):
        self.pixmap.draw_line(self.black_gc, int(x1), int(y1), int(x2), int(y2))






#
# DATA ------------------------------------------------------------------------------------------
#

    def data_set(self, data):
        self.data = data
        self.moviescale.set_sensitive(False)
        self.moviescale.set_range(0, 1)
        if len(self.data) > 1:
            self.moviescale.set_sensitive(True)
            self.moviescale.set_range(0, len(self.data) - 1)
            self.moviescale.connect("value-changed", lambda w: self.queue_draw())
        self.gfx_center_atoms()
        self.event_exposed()

#
# MAIN ------------------------------------------------------------------------------------------
#

if __name__ == "__main__":
    pid = os.fork()
    if pid:
        os._exit(0)
    import io
    import sys
    q = atomview()
    if len(sys.argv) > 1:
        data = io.loadposcars(sys.argv[1])
        if len(data) < 1:
            data = io.loadcons(sys.argv[1])
        q.data_set(data)
    gtk.main()
