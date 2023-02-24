import numpy as np

def _acceptable(matrix_shape,coord):
    x,y = coord
    if x>=0 and y>=0:
        if x<=matrix_shape[0]-1 and y<=matrix_shape[1]-1:
            return True
    return False

class Px_Contour():
    def __init__(self,row,col,fam=None,type="mask"):
        self.row = row
        self.col = col
        self.fam = fam
        self.type = type
        self.bounds = []
    def coord(self):
        return [self.row,self.col]
    def get_bordering_coord(self):
        bord_coords = [[self.row-1,self.col],
                       [self.row,self.col-1],
                       [self.row,self.col+1],
                       [self.row+1,self.col]]
        return bord_coords
    def get_bordering_pixels(self):
        bord_coords = self.get_bordering_coord()
        boord_pxs = [Px_Contour(col=bc[0],row=bc[1],fam=self.fam) for bc in bord_coords]
        return boord_pxs
    def __eq__(self, __o: object):
        return __o.row==self.row and __o.col==self.col

def get_neighboring_coords(matrix, coords):
    # find the neighboring pixels to the given ones 
    # and to which "family" of pixel it belongs to
    # divided into one for the pixels and one for the neighbors
    #neighbors = [] #the contours pixels
    #mpix_fam  = [] #family of the single pixel 
    #neigh_fam = [] #family of the contour pixel
    index_fam   = np.arange(len(matrix)*len(matrix[1])).reshape(np.shape(matrix))
    Px_Cnt_list= []
    for coord in coords:
        row, col = coord
        pxC = Px_Contour(row,col,fam=index_fam[row][col])
        if pxC in Px_Cnt_list:
            for pxCi in Px_Cnt_list:
                if pxC==pxCi:
                    pxC.fam = pxCi.fam
                    break
        else:
            Px_Cnt_list.append(pxC)
        nwpix = pxC.get_bordering_pixels()
        for nwP in nwpix:
            if _acceptable(np.shape(matrix),nwP.coord()):
                if nwP.coord() not in np.array(coords).tolist():
                    nwP.type="bound"
                if nwP in Px_Cnt_list and nwP.type=="mask":
                    for pxCi in Px_Cnt_list:
                        if nwP==pxCi:
                            fam_to_change=pxCi.fam
                            pxCi.fam = nwP.fam
                            break
                    for pxCi in Px_Cnt_list:
                        if pxCi.fam==fam_to_change:
                            pxCi.fam = nwP.fam
                else:
                    Px_Cnt_list.append(nwP)
    return Px_Cnt_list 

def interpolate_2D(data,pix,order=0):
    # data:2D matrix
    # pix: coordinate of the pixel to interpolate
    # order: order of interpolation. O=mean between neighboor pixels
    Pix_Cont_list = get_neighboring_coords(data,pix)
    fml_list      =  list(set([px.fam for px in Pix_Cont_list]))
    for fam in fml_list:
        valpx_fm_i  = []
        maskpx_fm_i =[]
        for px in Pix_Cont_list:
            if px.fam==fam:
                if px.type=="bound":
                    row,col = px.coord()
                    valpx_fm_i.append(data[row][col])
                else:
                    maskpx_fm_i.append(px.coord())
        if order==0:
            interp_val = np.nanmean(valpx_fm_i)
        else:
            #pragma: no cover
            raise
        for i,j in maskpx_fm_i:
            data[i][j] = interp_val
    return data 
