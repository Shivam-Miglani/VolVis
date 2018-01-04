/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */

// You do not need to change this class.
// This class represents the elements of the gradient volume.

public class VoxelGradient {

    public float x, y, z;
    public float mag;
    
    public VoxelGradient() {
        x = y = z = mag = 0.0f;
    }
    
    public VoxelGradient(float gx, float gy, float gz) {
        x = gx;
        y = gy;
        z = gz;
        mag = (float) Math.sqrt(x*x + y*y + z*z);
    }
    
}
