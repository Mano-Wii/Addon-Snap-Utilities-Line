### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 3
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ##### END GPL LICENSE BLOCK #####

# Contact for more information about the Addon:
# Email:    germano.costa@ig.com.br
# Twitter:  wii_mano @mano_wii

bl_info = {
    "name": "Snap_Utilities_Line",
    "author": "Germano Cavalcante",
    "version": (3, 0),
    "blender": (2, 74, 0),
    "location": "View3D > TOOLS > Snap Utilities > snap utilities",
    "description": "Extends Blender Snap controls",
    "wiki_url" : "http://blenderartists.org/forum/showthread.php?363859-Addon-CAD-Snap-Utilities",
    "category": "Mesh"}
    
import bpy, bgl, bmesh, mathutils, math
from mathutils import Vector, Matrix
from bpy_extras import view3d_utils

def location_3d_to_region_2d(region, rv3d, coord):
    prj = rv3d.perspective_matrix * Vector((coord[0], coord[1], coord[2], 1.0))
    width_half = region.width / 2.0
    height_half = region.height / 2.0
    return Vector((width_half + width_half * (prj.x / prj.w),
                   height_half + height_half * (prj.y / prj.w),
                   ))

def unProject(region, rv3d, mcursor):
    view_vector = view3d_utils.region_2d_to_vector_3d(region, rv3d, mcursor)
    ray_origin = view3d_utils.region_2d_to_origin_3d(region, rv3d, mcursor)
    ray_target = ray_origin + (view_vector * 1000)    
    scene = bpy.context.scene

    # cast the ray
    result, object, matrix, location, normal = scene.ray_cast(ray_origin, ray_target)
    if location == None:
        location = Vector((0,0,0))

    return result, object, matrix, location, normal

def out_Location(rv3d, region, mcursor):
    view_matrix = rv3d.view_matrix.transposed()
    orig = view3d_utils.region_2d_to_origin_3d(region, rv3d, mcursor)
    vector = view3d_utils.region_2d_to_vector_3d(region, rv3d, mcursor)
    v1 = Vector((int(view_matrix[0][0]*1.5),int(view_matrix[1][0]*1.5),int(view_matrix[2][0]*1.5)))
    v2 = Vector((int(view_matrix[0][1]*1.5),int(view_matrix[1][1]*1.5),int(view_matrix[2][1]*1.5)))
    
    hit = mathutils.geometry.intersect_ray_tri(Vector((1,0,0)), Vector((0,1,0)), Vector((0,0,0)), (vector), (orig), False)
    if hit == None:
        hit = mathutils.geometry.intersect_ray_tri(v1, v2, Vector((0,0,0)), (vector), (orig), False)        
    if hit == None:
        hit = mathutils.geometry.intersect_ray_tri(v1, v2, Vector((0,0,0)), (-vector), (orig), False)
    if hit == None:
        hit = Vector((0,0,0))
    return hit

def SnapUtilities(self, obj_matrix_world, bm_geom, bool_update, vert_perp, mcursor2, bool_constrain, vector_constrain):
    if not hasattr(self, 'const'):
        self.const = None

    if bool_constrain == False and self.const != None:
        self.const = None

    if isinstance(bm_geom, bmesh.types.BMVert):
        if not hasattr(self, 'bvert') or self.bvert != bm_geom or bool_update == True:
            self.bvert = bm_geom
            self.vert = obj_matrix_world * self.bvert.co
            self.Pvert = location_3d_to_region_2d(self.region, self.rv3d, self.vert)

        if bool_constrain == True:
            if self.const == None:
                self.const = self.vert
            #point = Vector([(self.vert[index] if vector_constrain==1 else self.const[index]) for index, vector_constrain in enumerate(vector_constrain)])
            point = mathutils.geometry.intersect_point_line(self.vert, self.const, (self.const+vector_constrain))[0]
            #point = vector_constrain.project(self.vert)
            return point, 'VERT' #50% is 'OUT'
        #else:
        return self.vert, 'VERT'
                
    if isinstance(bm_geom, bmesh.types.BMEdge):
        if not hasattr(self, 'bedge') or self.bedge != bm_geom or bool_update == True:
            self.bedge = bm_geom
            self.vert0 = obj_matrix_world*self.bedge.verts[0].co
            self.vert1 = obj_matrix_world*self.bedge.verts[1].co
            self.po_cent = (self.vert0+self.vert1)/2
            self.Pcent = location_3d_to_region_2d(self.region, self.rv3d, self.po_cent)
            self.Pvert0 = location_3d_to_region_2d(self.region, self.rv3d, self.vert0)
            self.Pvert1 = location_3d_to_region_2d(self.region, self.rv3d, self.vert1)
                
            if vert_perp != None and vert_perp not in [v.co for v in self.bedge.verts]:
                point_perpendicular = mathutils.geometry.intersect_point_line(vert_perp, self.vert0, self.vert1)
                self.po_perp = point_perpendicular[0]
                self.Pperp = location_3d_to_region_2d(self.region, self.rv3d, self.po_perp)

        if bool_constrain == True:
            if self.const == None:
                if vert_perp != None:
                    self.const = vert_perp
                else:
                    self.const = self.po_cent

            point = mathutils.geometry.intersect_line_line(self.const, (self.const+vector_constrain), self.vert0, self.vert1)
            if point == None:
                orig = view3d_utils.region_2d_to_origin_3d(self.region, self.rv3d, mcursor2 )
                view_vector = view3d_utils.region_2d_to_vector_3d(self.region, self.rv3d, mcursor2 )
                end = orig + view_vector
                point = mathutils.geometry.intersect_line_line(self.const, (self.const+vector_constrain), orig, end)
            return point[0], 'EDGE'

        else:
            if hasattr(self, 'Pperp') and abs(self.Pperp[0]-mcursor2[0]) < 10 and abs(self.Pperp[1]-mcursor2[1]) < 10:
                return self.po_perp, 'PERPENDICULAR'

            elif abs(self.Pcent[0]-mcursor2[0]) < 10 and abs(self.Pcent[1]-mcursor2[1]) < 10:
                return self.po_cent, 'CENTER'
            
            else:
                orig = view3d_utils.region_2d_to_origin_3d(self.region, self.rv3d, mcursor2)
                view_vector = view3d_utils.region_2d_to_vector_3d(self.region, self.rv3d, mcursor2)
                end = orig + view_vector
                point = mathutils.geometry.intersect_line_line(self.vert0, self.vert1, orig, end)
                return point[0], 'EDGE'

    if isinstance(bm_geom, bmesh.types.BMFace):
        if not hasattr(self, 'bface') or self.bface != bm_geom or bool_update == True:
            self.bface = bm_geom
            self.face_center = obj_matrix_world*bm_geom.calc_center_median()
            self.face_normal = bm_geom.normal*obj_matrix_world.inverted()
            
        orig = view3d_utils.region_2d_to_origin_3d(self.region, self.rv3d, mcursor2)
        view_vector = view3d_utils.region_2d_to_vector_3d(self.region, self.rv3d, mcursor2)
        end = orig + view_vector
        if bool_constrain == True:
            if self.const == None:
                if vert_perp != None:
                    self.const = vert_perp
                else:
                    self.const = point
            point = mathutils.geometry.intersect_line_line(self.const, (self.const+vector_constrain), orig, end)
            return point[0], 'FACE'
        #else:
        point = mathutils.geometry.intersect_line_plane(orig, end, self.face_center, self.face_normal, False)
        return point, 'FACE'
    
    else:
        if bool_constrain == True:
            if self.const == None:
                if vert_perp != None:
                    self.const = vert_perp
                else:
                    self.const = out_Location(self.rv3d, self.region, mcursor2)

            orig = view3d_utils.region_2d_to_origin_3d(self.region, self.rv3d, mcursor2 )
            view_vector = view3d_utils.region_2d_to_vector_3d(self.region, self.rv3d, mcursor2 )
            end = orig + view_vector
            point = mathutils.geometry.intersect_line_line(self.const, (self.const+vector_constrain), orig, end)
            return point[0], 'OUT'
        else:
            result, object, matrix, location, normal = unProject(self.region, self.rv3d, mcursor2)
            if result:
                return location, 'FACE'
            else:
                return out_Location(self.rv3d, self.region, mcursor2), 'OUT'

def get_isolated_edges(bmvert):
    linked = [c for c in bmvert.link_edges[:] if c.link_faces[:] == []]
    for a in linked:
        edges = [b for c in a.verts[:] if c.link_faces[:] == [] for b in c.link_edges[:] if b not in linked]
        for e in edges:
            linked.append(e)
    return linked

def split_face(self, mesh, Bmesh, listverts, listedges, listfaces):
    if len(listverts) >=2 and listverts[-1] not in [x for y in [a.verts[:] for a in listverts[-2].link_edges] for x in y if x != listverts[-2]]:
        if listverts[-2].link_faces[:] == []:
            if listfaces == []:
                for face in listverts[-1].link_faces:
                    testface = bmesh.geometry.intersect_face_point(face, listverts[-2].co)
                    if testface:
                        listfaces.append(face)
        else:
            face = [x for x in listverts[-1].link_faces[:] if x in listverts[-2].link_faces[:]]
            if face != []:
                listfaces.append(face[0])
        if listverts[-1] != listverts[-2]:
            facesp = {}
            facesp['edges'] = []
            if self.intersect and listverts[-2].link_faces[:] != []:
                verts = [listverts[-2], listverts[-1]]
                facesp = bmesh.ops.connect_vert_pair(Bmesh, verts = verts)
                bmesh.update_edit_mesh(mesh, tessface=True, destructive=True)
                #print(facesp)
            if not self.intersect or facesp['edges'] == []:
                edge = Bmesh.edges.new([listverts[-1], listverts[-2]])
                listedges.append(edge)
                if listfaces == []:
                    for face in listverts[-2].link_faces:
                        testface = bmesh.geometry.intersect_face_point(face, listverts[-1].co)
                        if testface:
                            listfaces.append(face)
                if listfaces != []:
                    for face in list(set(listfaces)):
                        facesp = bmesh.utils.face_split_edgenet(face, list(set(listedges)))
                        bmesh.update_edit_mesh(mesh, tessface=True, destructive=True)
            listedges = []

def draw_line(self, obj, Bmesh, bm_geom, location, bool_merge):
    if not hasattr(self, 'list_vertices'):
        self.list_vertices = []

    if not hasattr(self, 'list_edges'):
        self.list_edges = []

    if not hasattr(self, 'list_faces'):
        self.list_faces = []

    if bool_merge == False:
        vertices = (bmesh.ops.create_vert(Bmesh, co=(location)))
        self.list_vertices.append(vertices['vert'][0])
        split_face(self, obj.data, Bmesh, self.list_vertices, self.list_edges, self.list_faces)

    elif isinstance(bm_geom, bmesh.types.BMVert):
        if (bm_geom.co - location).length < .01:
            self.list_vertices.append(bm_geom)
            for edge in get_isolated_edges(bm_geom):
                if edge not in self.list_edges:
                    self.list_edges.append(edge)
        else:
            vertices = bmesh.ops.create_vert(Bmesh, co=(location))
            self.list_vertices.append(vertices['vert'][0])
            
        split_face(self, obj.data, Bmesh, self.list_vertices, self.list_edges, self.list_faces)
        
    elif isinstance(bm_geom, bmesh.types.BMEdge):
        self.list_edges.append(bm_geom)
        vector_p0_l = (bm_geom.verts[0].co-location)
        vector_p1_l = (bm_geom.verts[1].co-location)
        vector_p0_p1 = (bm_geom.verts[0].co-bm_geom.verts[1].co)
        factor = vector_p0_l.length/bm_geom.calc_length()
        if factor < 1 and vector_p1_l <  vector_p0_p1 and round(vector_p0_l.angle(vector_p1_l), 2) == 3.14: # contrain near                 
            vertex0 = bmesh.utils.edge_split(bm_geom, bm_geom.verts[0], factor)
            self.list_vertices.append(vertex0[1])
            self.list_edges.append(vertex0[0])
                                
            split_face(self, obj.data, Bmesh ,self.list_vertices, self.list_edges, self.list_faces)
            self.list_edges = []
        else:
            vertices = bmesh.ops.create_vert(Bmesh, co=(location))
            self.list_vertices.append(vertices['vert'][0])
            split_face(self, obj.data, Bmesh, self.list_vertices, self.list_edges, self.list_faces)

    elif isinstance(bm_geom, bmesh.types.BMFace):
        vertices = (bmesh.ops.create_vert(Bmesh, co=(location)))
        self.list_vertices.append(vertices['vert'][0])
        self.list_faces.append(bm_geom)
        split_face(self, obj.data, Bmesh, self.list_vertices, self.list_edges, self.list_faces)

    return [obj.matrix_world*a.co for a in self.list_vertices]

def draw_callback_px(self, context):
    # draw 3d point OpenGL in the 3D View
    bgl.glEnable(bgl.GL_BLEND)
    if self.bool_constrain:
        if self.vector_constrain == Vector((1,0,0)):
            Color4f = (self.axis_x_color + (1.0,))
        elif self.vector_constrain == Vector((0,1,0)):
            Color4f = (self.axis_y_color + (1.0,))
        elif self.vector_constrain == Vector((0,0,1)):
            Color4f = (self.axis_z_color + (1.0,))
        else:
            Color4f = self.constrain_shift_color
    else:
        if self.type == 'OUT':
            Color4f = self.out_color 
        elif self.type == 'FACE':
            Color4f = self.face_color
        elif self.type == 'EDGE':
            Color4f = self.edge_color
        elif self.type == 'VERT':
            Color4f = self.vert_color
        elif self.type == 'CENTER':
            Color4f = self.center_color
        elif self.type == 'PERPENDICULAR':
            Color4f = self.perpendicular_color

    bgl.glColor4f(*Color4f)
    bgl.glDepthRange(0,0)    
    bgl.glPointSize(10)    
    bgl.glBegin(bgl.GL_POINTS)
    bgl.glVertex3f(*self.location)
    bgl.glEnd()
    bgl.glDisable(bgl.GL_BLEND)

    # draw 3d line OpenGL in the 3D View
    bgl.glEnable(bgl.GL_BLEND)
    bgl.glDepthRange(0,0.9999)
    bgl.glColor4f(1.0, 0.8, 0.0, 1.0)    
    bgl.glLineWidth(2)    
    bgl.glEnable(bgl.GL_LINE_STIPPLE)
    bgl.glBegin(bgl.GL_LINE_STRIP)
    for vert_co in self.list_vertices_co:
        bgl.glVertex3f(*vert_co)        
    bgl.glVertex3f(*self.location)        
    bgl.glEnd()

    # restore opengl defaults
    bgl.glDepthRange(0,1)
    bgl.glPointSize(1)
    bgl.glLineWidth(1)
    bgl.glDisable(bgl.GL_BLEND)
    bgl.glDisable(bgl.GL_LINE_STIPPLE)
    bgl.glColor4f(0.0, 0.0, 0.0, 1.0)

    a = ""
    if self.list_vertices_co != [] and self.length_entered == "":
        a = 'length: '+ str(round((self.list_vertices_co[-1]-self.location).length, 3))
    elif self.list_vertices_co != [] and self.length_entered != "":
        a = 'length: '+ self.length_entered

    context.area.header_text_set("hit: %.3f %.3f %.3f %s" % (self.location[0], self.location[1], self.location[2], a))
    
class PanelSnapUtilities(bpy.types.Panel) :
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    #bl_context = "mesh_edit"
    bl_category = "Snap Utilities"
    bl_label = "snap utilities"
    
    @classmethod
    def poll(cls, context):
        return (context.object is not None and
                context.object.type == 'MESH')

    def draw(self, context):
        layout = self.layout
        TheCol = layout.column(align = True)
        TheCol.operator("mesh.snap_utilities_line", text = "Line", icon="GREASEPENCIL")
        
        addon_prefs = context.user_preferences.addons[__name__].preferences
        
        box = layout.box()
        row = box.row()
        row.prop(addon_prefs, "intersect")

class Constrain:
    keys = {
        'X': Vector((1,0,0)),
        'Y': Vector((0,1,0)),
        'Z': Vector((0,0,1)),
        'RIGHT_SHIFT': 'shift',
        'LEFT_SHIFT': 'shift',
        }

    def __init__(self, bool_constrain = False, vector_constrain = None):
        self.bool_constrain = bool_constrain
        self.vector_constrain = vector_constrain

    def modal(self, context, event):
        if event.value == 'PRESS':
            if self.vector_constrain == self.keys[event.type] or self.bool_constrain == False:
                self.bool_constrain = self.bool_constrain == False
                self.vector_constrain = self.keys[event.type]
                
            elif event.shift:
                if self.vector_constrain not in self.keys.values():
                    self.bool_constrain = self.bool_constrain == False
                    self.vector_constrain = self.keys[event.type]
                    
            else:
                self.vector_constrain = self.keys[event.type]
                    
        return self.bool_constrain, self.vector_constrain
    
class CharMap:
    keys = {
        'PERIOD':".", 'NUMPAD_PERIOD':".",
        'MINUS':"-", 'NUMPAD_MINUS':"-",
        'EQUAL':"+", 'NUMPAD_PLUS':"+",
        'ONE':"1", 'NUMPAD_1':"1",
        'TWO':"2", 'NUMPAD_2':"2",
        'THREE':"3", 'NUMPAD_3':"3",
        'FOUR':"4", 'NUMPAD_4':"4",
        'FIVE':"5", 'NUMPAD_5':"5",
        'SIX':"6", 'NUMPAD_6':"6",
        'SEVEN':"7", 'NUMPAD_7':"7",
        'EIGHT':"8", 'NUMPAD_8':"8",
        'NINE':"9", 'NUMPAD_9':"9",
        'ZERO':"0", 'NUMPAD_0':"0",
        'SPACE':" ",
        'SLASH':"/", 'NUMPAD_SLASH':"/",
        'NUMPAD_ASTERIX':"*",
        'BACK_SPACE':"", 'DEL':""
        }

    def __init__(self, length_entered = ""):
        self.length_entered = length_entered

    def modal(self, context, event):
        # Currently accessing event.ascii seems to crash Blender
        c = self.keys[event.type]
        if event.shift:
            if c == "8":
                c = "*"
            elif c == "5":
                c = "%"
            elif c == "9":
                c = "("
            elif c == "0":
                c = ")"

        if event.value == 'PRESS':
            self.length_entered += c
            if event.type in {'BACK_SPACE', 'DEL'} and len(self.length_entered) >= 1:
                self.length_entered = self.length_entered[:-1]

        return self.length_entered

class Navigation:
    keys = [
        'MIDDLEMOUSE',
        'WHEELDOWNMOUSE',
        'WHEELUPMOUSE',
        ]

    def __init__(self, rv3d, location):
        self.rv3d = rv3d
        self.location = location

    def modal(self, context, event):
        if event.type == 'MIDDLEMOUSE':
            self.rotMat = self.rv3d.view_matrix.copy()
            if event.value == 'PRESS':
                if event.shift:
                    if self.bool_constrain and (1 not in self.vector_constrain):
                        self.bool_constrain = False
                    bpy.ops.view3d.move('INVOKE_DEFAULT')
                        
                else:
                    bpy.ops.view3d.rotate('INVOKE_DEFAULT')

        if event.type in {'WHEELDOWNMOUSE', 'WHEELUPMOUSE'}:
            delta = (event.type == 'WHEELUPMOUSE') - (event.type == 'WHEELDOWNMOUSE')
            self.rv3d.view_distance -= delta*self.rv3d.view_distance/6
            self.rv3d.view_location += delta*(self.location - self.rv3d.view_location)/6

class MESH_OT_snap_utilities_line(bpy.types.Operator):
    """ Draw edges. Connect them to split faces."""
    bl_idname = "mesh.snap_utilities_line"
    bl_label = "Line Tool"
    bl_options = {'REGISTER', 'UNDO'}
    
    def modal(self, context, event):
        if context.area:
            context.area.tag_redraw()
            if self.rv3d.view_matrix != self.rotMat:
                self.rotMat = self.rv3d.view_matrix
                self.bool_update = True
            else:
                self.bool_update = False
                
        if event.ctrl:
            if event.type == 'Z' and event.value == 'PRESS':
                bpy.ops.ed.undo()
                self.bool_constrain = False
                self.list_vertices_co = []
                self.list_vertices = []
                self.list_edges = []
                self.list_faces = []
                self.obj = bpy.context.active_object
                self.obj_matrix = self.obj.matrix_world.copy()
                self.bm = bmesh.from_edit_mesh(self.obj.data)
                return {'RUNNING_MODAL'}

        if event.type in Navigation.keys:
            Navigation.modal(self, context, event)

        elif event.type in Constrain.keys:
            Constrain2 = Constrain(self.bool_constrain, self.vector_constrain)
            self.bool_constrain, self.vector_constrain = Constrain2.modal(context, event)
            if self.vector_constrain == 'shift':
                if isinstance(self.geom, bmesh.types.BMEdge):
                    #self.vector_constrain = self.obj_matrix*self.geom.verts[1].co-self.obj_matrix*self.geom.verts[0].co
                    self.vector_constrain = (self.geom.verts[1].co-self.geom.verts[0].co)*self.obj_matrix.inverted()
                else:
                    self.bool_constrain = False

        elif event.type in CharMap.keys and event.value == 'PRESS':
            CharMap2 = CharMap(self.length_entered)
            self.length_entered = CharMap2.modal(context, event)
            #print(self.length_entered)
        
        elif event.type == 'MOUSEMOVE':
            x, y = (event.mouse_region_x, event.mouse_region_y)
            if self.obj:
                bpy.ops.mesh.select_all(action='DESELECT')
                
            bpy.ops.view3d.select(location=(x, y))

            if self.list_vertices_co != []:
                bm_vert_to_perpendicular = self.list_vertices_co[-1]
            else:
                bm_vert_to_perpendicular = None
            
            try:
                self.geom = self.bm.select_history[0]
            except: # IndexError or AttributeError:
                self.geom = None

            self.location, self.type = SnapUtilities(self, self.obj_matrix, self.geom, self.bool_update, bm_vert_to_perpendicular, (x, y), self.bool_constrain, self.vector_constrain)
            
        elif event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            #if event.value == 'PRESS':
            # SNAP 2D
            snap_3d = self.location
            Lsnap_3d = self.obj_matrix.inverted()*snap_3d
            Snap_2d = location_3d_to_region_2d(self.region, self.rv3d, snap_3d)
            if self.bool_constrain and isinstance(self.geom, bmesh.types.BMVert): # SELECT FIRST
                bpy.ops.view3d.select(location=(int(Snap_2d[0]), int(Snap_2d[1])))
                try:
                    geom2 = self.bm.select_history[0]
                except: # IndexError or AttributeError:
                    geom2 = None
            else:
                geom2 = self.geom
            bool_merge = self.type not in {'FACE', 'OUT'}
            self.bool_constrain = False
            self.list_vertices_co = draw_line(self, self.obj, self.bm, geom2, Lsnap_3d, bool_merge)
            bpy.ops.ed.undo_push(message="Add an undo step *function may be moved*")
            
        elif event.type in {'RET', 'NUMPAD_ENTER'} and event.value == 'RELEASE':
            if self.length_entered != "" and self.list_vertices_co != []:
                try:
                    text_value = eval(self.length_entered, math.__dict__)
                    vector_h0_h1 = (self.location-self.list_vertices_co[-1]).normalized()
                    location = ((vector_h0_h1*text_value)+self.obj_matrix.inverted()*self.list_vertices_co[-1])
                    bool_merge = self.type not in {'FACE', 'OUT'}
                    self.list_vertices_co = draw_line(self, self.obj, self.bm, self.geom, location, bool_merge)
                    self.length_entered = ""
                
                except:# ValueError:
                    self.report({'INFO'}, "Operation not supported yet")
        
        elif event.type == 'TAB' and event.value == 'PRESS':
            self.keytab = self.keytab == False
            if self.keytab:            
                context.tool_settings.mesh_select_mode = (False, False, True)
            else:
                context.tool_settings.mesh_select_mode = (True, True, True)
            
        elif event.type in {'RIGHTMOUSE', 'ESC'} and event.value == 'RELEASE':
            if self.list_vertices_co == [] or event.type == 'ESC':                
                bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
                context.tool_settings.mesh_select_mode = self.select_mode
                context.area.header_text_set()
                context.user_preferences.view.use_rotate_around_active = self.use_rotate_around_active
                if not self.is_editmode:
                    bpy.ops.object.editmode_toggle()
                return {'FINISHED'}
            else:
                self.bool_constrain = False
                self.list_vertices = []
                self.list_vertices_co = []
                self.list_faces = []
                
        return {'RUNNING_MODAL'}

    def invoke(self, context, event):        
        if context.space_data.type == 'VIEW_3D':
            bgl.glEnable(bgl.GL_POINT_SMOOTH)
            self.is_editmode = bpy.context.object.data.is_editmode
            bpy.ops.object.mode_set(mode='EDIT')
            context.space_data.use_occlude_geometry = True

            self.use_rotate_around_active = context.user_preferences.view.use_rotate_around_active
            context.user_preferences.view.use_rotate_around_active = True
            
            self.select_mode = context.tool_settings.mesh_select_mode[:]
            context.tool_settings.mesh_select_mode = (True, True, True)
            
            self.region = context.region
            self.rv3d = context.region_data
            self.rotMat = self.rv3d.view_matrix
            self.obj = bpy.context.active_object
            self.obj_matrix = self.obj.matrix_world.copy()
            self.bm = bmesh.from_edit_mesh(self.obj.data)
            
            self.list_vertices_co = []
            self.bool_constrain = False
            self.bool_update = False
            self.vector_constrain = None
            self.keytab = False
            self.length_entered = ""
            self._handle = bpy.types.SpaceView3D.draw_handler_add(draw_callback_px, (self, context), 'WINDOW', 'POST_VIEW')
            context.window_manager.modal_handler_add(self)
            
            self.out_color = context.user_preferences.addons[__name__].preferences.out_color
            self.face_color = context.user_preferences.addons[__name__].preferences.face_color
            self.edge_color = context.user_preferences.addons[__name__].preferences.edge_color
            self.vert_color = context.user_preferences.addons[__name__].preferences.vert_color
            self.center_color = context.user_preferences.addons[__name__].preferences.center_color
            self.perpendicular_color = context.user_preferences.addons[__name__].preferences.perpendicular_color
            self.constrain_shift_color = context.user_preferences.addons[__name__].preferences.constrain_shift_color

            self.axis_x_color = tuple(context.user_preferences.themes[0].user_interface.axis_x)
            self.axis_y_color = tuple(context.user_preferences.themes[0].user_interface.axis_y)
            self.axis_z_color = tuple(context.user_preferences.themes[0].user_interface.axis_z)

            self.intersect = context.user_preferences.addons[__name__].preferences.intersect

            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}

class SnapAddonPreferences(bpy.types.AddonPreferences):
    # this must match the addon name, use '__package__'
    # when defining this in a submodule of a python package.
    bl_idname = __name__
    
    intersect = bpy.props.BoolProperty(
            name="intersect",
            description="intersects created line with the existing edges, even if the lines do not intersect.",
            default=True,
            )
    out_color = bpy.props.FloatVectorProperty(name="OUT", default=(0.0, 0.0, 0.0, 0.5), size=4, subtype="COLOR", min=0, max=1)
    face_color = bpy.props.FloatVectorProperty(name="FACE", default=(1.0, 0.8, 0.0, 1.0), size=4, subtype="COLOR", min=0, max=1)
    edge_color = bpy.props.FloatVectorProperty(name="EDGE", default=(0.0, 0.8, 1.0, 1.0), size=4, subtype="COLOR", min=0, max=1)
    vert_color = bpy.props.FloatVectorProperty(name="VERT", default=(1.0, 0.5, 0.0, 1.0), size=4, subtype="COLOR", min=0, max=1)
    center_color = bpy.props.FloatVectorProperty(name="CENTER", default=(1.0, 0.0, 1.0, 1.0), size=4, subtype="COLOR", min=0, max=1)
    perpendicular_color = bpy.props.FloatVectorProperty(name="PERPENDICULAR", default=(0.1, 0.5, 0.5, 1.0), size=4, subtype="COLOR", min=0, max=1)
    constrain_shift_color = bpy.props.FloatVectorProperty(name="SHIFT CONSTRAIN", default=(0.8, 0.5, 0.4, 1.0), size=4, subtype="COLOR", min=0, max=1)

    def draw(self, context):
        layout = self.layout

        layout.label(text="Snap Colors:")
        split = layout.split()

        col = split.column()
        col.prop(self, "out_color")
        col.prop(self, "constrain_shift_color")
        col = split.column()
        col.prop(self, "face_color")
        col = split.column()
        col.prop(self, "edge_color")        
        col = split.column()
        col.prop(self, "vert_color")
        col = split.column()
        col.prop(self, "center_color")
        col = split.column()
        col.prop(self, "perpendicular_color")
        
        layout.label(text="Line tool:")        
        layout.prop(self, "intersect")

def register():
    print('Addon', __name__, 'registered')
    bpy.utils.register_class(SnapAddonPreferences)
    bpy.utils.register_class(PanelSnapUtilities)
    bpy.utils.register_class(MESH_OT_snap_utilities_line)

def unregister():
    bpy.utils.unregister_class(SnapAddonPreferences)
    bpy.utils.unregister_class(PanelSnapUtilities)
    bpy.utils.unregister_class(MESH_OT_snap_utilities_line)

if __name__ == "__main__":
    __name__ = "mesh_snap_utilities_line"
    register()
