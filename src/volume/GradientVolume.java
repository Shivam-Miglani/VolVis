/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 *  Modified by Anna Vilanova
 */
public class GradientVolume {

	
	
	//////////////////////////////////////////////////////////////////////
	///////////////// TO BE IMPLEMENTED //////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
	//Compute the gradient of contained in the volume attribute and save it into the data attribute
	//
	//This is a lengthy computation and is performed only once (have a look at the constructor GradientVolume) 

    // You need to implement this function
    private void compute() {
        int i, j, k;
        float gx, gy, gz;
        for (int q =0; q< data.length; q++){
            data[q]= new VoxelGradient();
        }
        
        
        for(i=1;i<dimX-1;i++){
            for(j=1;j<dimY-1;j++){
                for(k=1;k<dimZ-1;k++){
                    
                    gx= (volume.data[(i+1)+(dimX*(j+(dimY*k)))] - volume.data[(i-1)+(dimX*(j+(dimY*k)))]) / 2;
                    
                    gy= (volume.data[i+(dimX*((j+1)+(dimY*k)))] - volume.data[i+(dimX*((j-1)+(dimY*k)))]) / 2;
                    
                    gz= (volume.data[i+(dimX*(j+(dimY*(k+1))))] - volume.data[i+(dimX*(j+(dimY*(k-1))))]) / 2;
                    
                    data[i + dimX * (j + dimY * k)]= new VoxelGradient(gx, gy, gz);
                    
                    
                }
            }
        }

        // to be implemented
    }
    	
    //You need to implement this function
    //This function linearly interpolates gradient vector g0 and g1 given the factor (t) 
    //the resut is given at result. You can use it to tri-linearly interpolate the gradient
    private VoxelGradient interpolate(VoxelGradient g0, VoxelGradient g1, float factor) {
       // to be implemented
       float x, y, z;
       VoxelGradient result;
       
       x= (g0.x*(1-factor)) + (g1.x*factor);
       y= (g0.y*(1-factor)) + (g1.y*factor);
       z= (g0.z*(1-factor)) + (g1.z*factor);
       
       result= new VoxelGradient(x, y, z);
       
       return result;
    }
    
    // You need to implement this function
    // This function returns the gradient at position coord using trilinear interpolation
    public VoxelGradient getGradient(double[] coord) {

        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }
        
        float x = (float)coord[0]; 
        float y = (float)coord[1];
        float z = (float)coord[2];
        
        float x0 = (float) Math.floor(x);
        float y0 = (float) Math.floor(y);
        float z0 = (float) Math.floor(z);
        float x1 = (float) Math.ceil(x);
        float y1 = (float) Math.ceil(y);
        float z1 = (float) Math.ceil(z);
        
        float xfactor = (x-x0)/(x1-x0);
        float yfactor = (y-y0)/(y1-y0);
        float zfactor = (z-z0)/(z1-z0);
        
        VoxelGradient c000 = new VoxelGradient (x0,y0,z0);
        VoxelGradient c001 = new VoxelGradient (x0,y0,z1);
        VoxelGradient c010 = new VoxelGradient (x0,y1,z0);
        VoxelGradient c011 = new VoxelGradient (x0,y1,z1);
        VoxelGradient c100 = new VoxelGradient (x1,y0,z0);
        VoxelGradient c101 = new VoxelGradient (x1,y0,z1);
        VoxelGradient c110 = new VoxelGradient (x1,y1,z0);
        VoxelGradient c111 = new VoxelGradient (x1,y1,z1);
        
        VoxelGradient c00= interpolate(c000,c100,xfactor);
        VoxelGradient c01= interpolate(c001,c101,xfactor);
        VoxelGradient c10= interpolate(c010,c110,xfactor);
        VoxelGradient c11= interpolate(c011,c111,xfactor);
        
        VoxelGradient c0= interpolate(c00,c10,yfactor);
        VoxelGradient c1= interpolate(c01,c11,yfactor);
        
        VoxelGradient c = interpolate(c0,c1,zfactor);
        
        //x= c.x;
        //y= c.y;
        //z= c.z;

        //int x = (int) Math.floor(coord[0]);
        //int y = (int) Math.floor(coord[1]);
        //int z = (int) Math.floor(coord[2]);
        //System.out.println(getGradient((int)x, (int)y, (int)z).mag);
        //return getGradient((int)x, (int)y, (int)z);
        return c;
        // To be impmented right now it impments just a nearest neighbour
    }
    
    
    
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
	
    //Returns the maximum gradient magnitude
    //
    //The data array contains all the gradients, in this function you have to return the maximum magnitude of the vectors in data[] 
 
    //Do NOT modify this function
    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    	
	
	
	//Do NOT modify this function
	public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }

	//Do NOT modify this function
	public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }

  
	//Do NOT modify this function: Basically calculates the Nearest Neighbor interpolation for the gradient
    public VoxelGradient getGradient2(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }

        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
        return getGradient(x, y, z);
    }

	//Do NOT modify this function
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

	//Do NOT modify this function
    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }
    
	//Do NOT modify this function
    public VoxelGradient getVoxel(int i) {
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

	//Do NOT modify this attributes
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
    
    //If needed add new attributes here:
}
