/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

/**
 *
 * @author anna
 */
public class TransferFunction2D {
    public short baseIntensity;
    public double radius;
    public TFColor color;
        

        public TransferFunction2D(short base, double r) {
            this.baseIntensity = base;
            this.radius = r;
            this.color = new TFColor(0.0, 204.0/255.0, 153.0/255.0, 0.3);
        }
        
        public void SetBaseRadius(short base, double r)
        {   
            this.baseIntensity = base;
            this.radius = r;
        }
}
