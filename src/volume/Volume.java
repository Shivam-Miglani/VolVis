/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import java.io.File;
import java.io.IOException;

/**
 * @author michel
 *         Modified by Anna Vilanova
 *         Modified by Galib, Sharad and Shivam. (group 7)
 */
public class Volume {

    //This function linearly interpolates the value g0 and g1 given the factor (t) 
    //the result is returned. You can use it to tri-linearly interpolate the values 
    private float interpolate(float g0, float g1, float factor) {
        float result;
        // to be implemented
        result = (g0 * (1 - factor)) + (g1 * factor);

        return result;
    }

    private float[] vectorInterpolate(float[] c1, float[] c2, float factor) {
        float[] cnew = new float[3];
        for (int i = 0; i < c1.length; i++) {
            cnew[i] = interpolate(c1[i], c2[i], factor);
        }
        return cnew;
    }

    //You have to implement the trilinear interpolation of the volume
    //First implement the interpolated function above
    // At the moment the function does takes just the lowest voxel value
    // to trilinear interpolation
    public short getVoxelLinearInterpolate(double[] coord) {
        /* if coord are not in proper format, return 0*/
        if (coord[0] < 0 || coord[0] > (dimX - 2) || coord[1] < 0 || coord[1] > (dimY - 2)
                || coord[2] < 0 || coord[2] > (dimZ - 2)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        float x = (float) coord[0];
        float y = (float) coord[1];
        float z = (float) coord[2];

        int x0 = (int) Math.floor(x);
        int y0 = (int) Math.floor(y);
        int z0 = (int) Math.floor(z);
        int x1 = (int) Math.ceil(x);
        int y1 = (int) Math.ceil(y);
        int z1 = (int) Math.ceil(z);

        float xfactor = (x - x0) / (x1 - x0);
        float yfactor = (y - y0) / (y1 - y0);
        float zfactor = (z - z0) / (z1 - z0);

        /*
        *** The 7 interpolations ***
        */

        /* 4 linear interpolations in x direction*/
        float c00 = interpolate(data[x0 + dimX * (y0 + dimY * z0)], data[(x1) + dimX * (y0 + dimY * z0)], xfactor);
        float c01 = interpolate(data[x0 + dimX * (y0 + dimY * z1)], data[(x1) + dimX * (y0 + dimY * z1)], xfactor);
        float c10 = interpolate(data[x0 + dimX * (y1 + dimY * z0)], data[(x1) + dimX * (y1 + dimY * z0)], xfactor);
        float c11 = interpolate(data[x0 + dimX * (y1 + dimY * z1)], data[(x1) + dimX * (y1 + dimY * z1)], xfactor);

        /*2 linear interpolations in y dimension*/
        float c0 = interpolate(c00, c10, yfactor);
        float c1 = interpolate(c01, c11, yfactor);

        /*final linear interpolation in z direction*/
        float c = interpolate(c0, c1, zfactor);

        //System.out.println((short)c);
        return (short) c;
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    //Do NOT modify this function
    // This function is an example and does a nearest neighbour interpolation
    public short getVoxelNN(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX - 1) || coord[1] < 0 || coord[1] > (dimY - 1)
                || coord[2] < 0 || coord[2] > (dimZ - 1)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);

        return getVoxel(x, y, z);
    }

    //Do NOT modify this function
    public Volume(int xd, int yd, int zd) {
        data = new short[xd * yd * zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }

    //Do NOT modify this function
    public Volume(File file) {

        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            System.out.println(dimX + " " + dimY + " " + dimZ);
            System.out.println(data.length);
            System.out.println(data[35252]);
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }

    }

    //Do NOT modify this function
    public short getVoxel(int x, int y, int z) {
        int i = x + dimX * (y + dimY * z);
        return data[i];
    }

    //Do NOT modify this function
    public void setVoxel(int x, int y, int z, short value) {
        int i = x + dimX * (y + dimY * z);
        data[i] = value;
    }

    //Do NOT modify this function
    public void setVoxel(int i, short value) {
        data[i] = value;
    }

    //Do NOT modify this function
    public short getVoxel(int i) {
        return data[i];
    }

    //Do NOT modify this function
    public int getDimX() {
        return dimX;
    }

    //Do NOT modify this function
    public int getDimY() {
        return dimY;
    }

    //Do NOT modify this function
    public int getDimZ() {
        return dimZ;
    }

    //Do NOT modify this function
    public short getMinimum() {
        short minimum = data[0];
        for (int i = 0; i < data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }

    //Do NOT modify this function
    public short getMaximum() {
        short maximum = data[0];
        for (int i = 0; i < data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }

    //Do NOT modify this function
    public int[] getHistogram() {
        return histogram;
    }

    //Do NOT modify this function
    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i = 0; i < data.length; i++) {
            histogram[data[i]]++;
        }
    }

    //Do NOT modify these attributes
    private int dimX, dimY, dimZ;
    short[] data;
    private int[] histogram;
}
