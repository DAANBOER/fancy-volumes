/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import java.io.File;

import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

import javax.imageio.ImageIO;
import javax.xml.bind.util.ValidationEventCollector;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;

    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());

        // uncomment this to initialize the TF with good starting values for the orange dataset
        //tFunc.setTestFunc();


        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());

        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }

    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }


    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] >= volume.getDimX() || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }


    short getVoxelTriLinear(double[] coord) {
        if (coord[0] < 0 || coord[0] >= volume.getDimX() || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            return 0;
        }

        double[] alpha = {0, 0, 0};

        int[][] cubeCoords = {{0, 0}, {0, 0}, {0, 0}};

        for (int i = 0; i < 3; i++) {
            cubeCoords[i][0] = (int) Math.floor(coord[i]);
            cubeCoords[i][1] = (int) Math.ceil(coord[i]);
        }

        if (cubeCoords[0][1] >= volume.getDimX() || cubeCoords[1][1] >= volume.getDimY() || cubeCoords[2][1] >= volume.getDimZ()) {
            return 0;
        }

        for (int i = 0; i < 3; i++) {
            alpha[i] = (coord[i] - cubeCoords[i][0]) / (double)(cubeCoords[i][1] - cubeCoords[i][0]);
        }

        short[][] interX = {{0, 0}, {0, 0}};
        short[] interY = {0, 0};

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                interX[i][j] = (short) ( (1 - alpha[0]) * volume.getVoxel(cubeCoords[0][0], cubeCoords[1][i], cubeCoords[2][j]) +
                        alpha[0] * volume.getVoxel(cubeCoords[0][1], cubeCoords[1][i], cubeCoords[2][j]));
            }
        }

        interY[0] = (short) ((1 - alpha[1]) * interX[0][0] + alpha[1] * interX[1][0]);
        interY[1] = (short) ((1 - alpha[1]) * interX[0][1] + alpha[1] * interX[1][1]);

        return (short) ((1 - alpha[2]) * interY[0] + alpha[2] * interY[1]);
    }

    float getGradientMagnitudeTrilinear(double[] coord) {
        if (coord[0] < 0 || coord[0] >= gradients.getDimX() || coord[1] < 0 || coord[1] >= gradients.getDimY()
                || coord[2] < 0 || coord[2] >= gradients.getDimZ()) {
            return 0;
        }

        double[] alpha = {0, 0, 0};

        int[][] cubeCoords = {{0, 0}, {0, 0}, {0, 0}};

        for (int i = 0; i < 3; i++) {
            cubeCoords[i][0] = (int) Math.floor(coord[i]);
            cubeCoords[i][1] = (int) Math.ceil(coord[i]);
        }

        if (cubeCoords[0][1] >= gradients.getDimX() || cubeCoords[1][1] >= gradients.getDimY() || cubeCoords[2][1] >= gradients.getDimZ()) {
            return 0;
        }

        for (int i = 0; i < 3; i++) {
            alpha[i] = (coord[i] - cubeCoords[i][0]) / (double)(cubeCoords[i][1] - cubeCoords[i][0]);
        }

        float[][] interX = {{0, 0}, {0, 0}};
        float[] interY = {0, 0};

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                interX[i][j] = (float) ( (1. - alpha[0]) * gradients.getGradient(cubeCoords[0][0], cubeCoords[1][i], cubeCoords[2][j]).mag +
                        alpha[0] * gradients.getGradient(cubeCoords[0][1], cubeCoords[1][i], cubeCoords[2][j]).mag);
            }
        }

        interY[0] = (float) ((1. - alpha[1]) * interX[0][0] + alpha[1] * interX[1][0]);
        interY[1] = (float) ((1. - alpha[1]) * interX[0][1] + alpha[1] * interX[1][1]);

        return (float) ((1. - alpha[2]) * interY[0] + alpha[2] * interY[1]);
    }

    VoxelGradient getGradientTrilinear(double[] coord) {
        if (coord[0] < 0 || coord[0] >= gradients.getDimX() || coord[1] < 0 || coord[1] >= gradients.getDimY()
                || coord[2] < 0 || coord[2] >= gradients.getDimZ()) {
            return new VoxelGradient();
        }

        double[] alpha = {0, 0, 0};

        int[][] cubeCoords = {{0, 0}, {0, 0}, {0, 0}};

        for (int i = 0; i < 3; i++) {
            cubeCoords[i][0] = (int) Math.floor(coord[i]);
            cubeCoords[i][1] = (int) Math.ceil(coord[i]);
        }

        if (cubeCoords[0][1] >= gradients.getDimX() || cubeCoords[1][1] >= gradients.getDimY() || cubeCoords[2][1] >= gradients.getDimZ()) {
            return new VoxelGradient();
        }

        for (int i = 0; i < 3; i++) {
            alpha[i] = (coord[i] - cubeCoords[i][0]) / (double)(cubeCoords[i][1] - cubeCoords[i][0]);
        }

        VoxelGradient[][] xGradients = new VoxelGradient[2][2];
        VoxelGradient[] yGradients = new VoxelGradient[2];

        float x, y, z;

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                x = (float) ( (1. - alpha[0]) * gradients.getGradient(cubeCoords[0][0], cubeCoords[1][i], cubeCoords[2][j]).x +
                        alpha[0] * gradients.getGradient(cubeCoords[0][1], cubeCoords[1][i], cubeCoords[2][j]).x);
                y = (float) ( (1. - alpha[0]) * gradients.getGradient(cubeCoords[0][0], cubeCoords[1][i], cubeCoords[2][j]).y +
                        alpha[0] * gradients.getGradient(cubeCoords[0][1], cubeCoords[1][i], cubeCoords[2][j]).y);
                z = (float) ( (1. - alpha[0]) * gradients.getGradient(cubeCoords[0][0], cubeCoords[1][i], cubeCoords[2][j]).z +
                        alpha[0] * gradients.getGradient(cubeCoords[0][1], cubeCoords[1][i], cubeCoords[2][j]).z);

                xGradients[i][j] = new VoxelGradient(x, y, z);
            }
        }

        for (int i = 0; i < 2; i++) {
            x = (float)((1. - alpha[1]) * xGradients[0][i].x + alpha[1] + xGradients[1][i].x);
            y = (float)((1. - alpha[1]) * xGradients[0][i].y + alpha[1] + xGradients[1][i].y);
            z = (float)((1. - alpha[1]) * xGradients[0][i].z + alpha[1] + xGradients[1][i].z);

            yGradients[i] = new VoxelGradient(x, y, z);
        }

        return new VoxelGradient((float) ((1. - alpha[2]) * yGradients[0].x + alpha[2] * yGradients[1].x),
                (float) ((1. - alpha[2]) * yGradients[0].y + alpha[2] * yGradients[1].y),
                (float) ((1. - alpha[2]) * yGradients[0].z + alpha[2] * yGradients[1].z));

        //return (float) ((1. - alpha[2]) * interY[0] + alpha[2] * interY[1]);
    }

    VoxelGradient getGradient(double[] coord) {
        if (coord[0] < 0 || coord[0] >= volume.getDimX() || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            return new VoxelGradient();
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return gradients.getGradient(x, y, z);
    }

    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, (double)volume.getDimX() / 2., (double)volume.getDimY() / 2., (double)volume.getDimZ() / 2.);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                //int val = getVoxel(pixelCoord);
                short val = getVoxelTriLinear(pixelCoord);

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);


                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }

    void colorSlicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, (double)volume.getDimX() / 2., (double)volume.getDimY() / 2., (double)volume.getDimZ() / 2.);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                //int val = getVoxel(pixelCoord);
                short val = getVoxelTriLinear(pixelCoord);

                voxelColor = tFunc.getColor(val);

                // Map the intensity to a grey value by linear scaling
                /*voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0; */ // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);


                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                image.setRGB(i, j, pixelColor);
            }
        }

    }

    void mip(double[] viewMatrix, int resolutionScaling, int depthScaling) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, (double)volume.getDimX() / 2., (double)volume.getDimY() / 2., (double)volume.getDimZ() / 2.);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        int viewDim = (int) Math.sqrt(Math.pow((double) volume.getDimZ(), 2.) +
                                      Math.pow((double) volume.getDimX(), 2.) +
                                      Math.pow((double) volume.getDimY(), 2.) );

        short val = 0, maxVal;

        for (int j = 0; j < image.getHeight(); j+=resolutionScaling) {
            for (int i = 0; i < image.getWidth(); i+=resolutionScaling) {

                /*viewDim = (int) Math.sqrt(Math.pow((double) volume.getDimZ(), 2.) +
                                          Math.pow((double) volume.getDimX(), 2.) +
                                          Math.pow((double) volume.getDimY(), 2.) );*/

                maxVal = Short.MIN_VALUE;

                for (int k = 0; k < viewDim; k+=depthScaling) {

                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];

                    try {
                        val = getVoxelTriLinear(pixelCoord);
                    } catch (ArrayIndexOutOfBoundsException e) {
                        System.out.println("trilinear is the problem");
                    }


                    if (val > maxVal) maxVal = val;
                }

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxVal/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxVal > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);


                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                setPixels(i, j, resolutionScaling, depthScaling, pixelColor);
            }
        }
    }

    void composite(double[] viewMatrix, int resolutionScaling, int depthScaling) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, (double)volume.getDimX() / 2., (double)volume.getDimY() / 2., (double)volume.getDimZ() / 2.);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        int viewDim = (int) Math.sqrt(Math.pow((double) volume.getDimZ(), 2.) +
                Math.pow((double) volume.getDimX(), 2.) +
                Math.pow((double) volume.getDimY(), 2.) );

        TFColor current, composite;

        for (int j = 0; j < image.getHeight(); j+=resolutionScaling) {
            for (int i = 0; i < image.getWidth(); i+=resolutionScaling) {

                //short val = 0;
                composite = new TFColor(0., 0., 0., 0.);

                for (int k = viewDim - 1; k >= 0; k-=depthScaling) {

                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];

                    current = tFunc.getColor(getVoxelTriLinear(pixelCoord));

                    composite.r = current.r * current.a + (1. - current.a) * composite.r;
                    composite.g = current.g * current.a + (1. - current.a) * composite.g;
                    composite.b = current.b * current.a + (1. - current.a) * composite.b;
                }

                // Map the intensity to a grey value by linear scaling

                voxelColor.r = composite.r;
                voxelColor.g = composite.g;
                voxelColor.b = composite.b;
                voxelColor.a = 1.;

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                setPixels(i, j, resolutionScaling, depthScaling, pixelColor);
            }
        }
    }

    void compositeFront(double[] viewMatrix, int resolutionScaling, int depthScaling) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, (double)volume.getDimX() / 2., (double)volume.getDimY() / 2., (double)volume.getDimZ() / 2.);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        int viewDim = (int) Math.sqrt(Math.pow((double) volume.getDimZ(), 2.) +
                Math.pow((double) volume.getDimX(), 2.) +
                Math.pow((double) volume.getDimY(), 2.) );

        TFColor current, composite;

        for (int j = 0; j < image.getHeight(); j+=resolutionScaling) {
            for (int i = 0; i < image.getWidth(); i+=resolutionScaling) {

                composite = new TFColor(0., 0., 0., 1.);

                for (int k = 0; k < viewDim; k += depthScaling) {

                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];

                    current = tFunc.getColor(getVoxelTriLinear(pixelCoord));

                    composite.r = composite.r + composite.a * current.a * current.r;
                    composite.g = composite.g + composite.a * current.a * current.g;
                    composite.b = composite.b + composite.a * current.a * current.b;
                    composite.a *= 1. - current.a;

                    if (composite.a < 0.05) {
                        break;
                    }
                }

                voxelColor.r = composite.r;
                voxelColor.g = composite.g;
                voxelColor.b = composite.b;
                voxelColor.a = 1. - composite.a;

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                setPixels(i, j, resolutionScaling, depthScaling, pixelColor);
            }
        }
    }

    void transfer2D(double[] viewMatrix, int resolutionScaling, int depthScaling) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, (double)volume.getDimX() / 2., (double)volume.getDimY() / 2., (double)volume.getDimZ() / 2.);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        int viewDim = (int) Math.sqrt(Math.pow((double) volume.getDimZ(), 2.) +
                Math.pow((double) volume.getDimX(), 2.) +
                Math.pow((double) volume.getDimY(), 2.) );

        short voxel;
        float gradient;
        double alpha, totalAlpha;

        double r = tfEditor2D.triangleWidget.radius;
        double av = tfEditor2D.triangleWidget.color.a;
        short fv = tfEditor2D.triangleWidget.baseIntensity;

        for (int j = 0; j < image.getHeight(); j+=resolutionScaling) {
            for (int i = 0; i < image.getWidth(); i+=resolutionScaling) {

                totalAlpha = 1.;

                for (int k = 0; k < viewDim; k += depthScaling) {

                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];

                    voxel = getVoxelTriLinear(pixelCoord);

                    gradient = getGradientMagnitudeTrilinear(pixelCoord);

                    if (gradient < tfEditor2D.triangleWidget.minGradMag || gradient > tfEditor2D.triangleWidget.maxGradMag) {
                        alpha = 0;
                    } else {
                        if (gradient == 0. && fv == voxel) {
                            alpha = av;
                        } else if (gradient > 0 && fv >= voxel - r * gradient && fv <= voxel + r * gradient) {
                            alpha = av * (1. - Math.abs((fv - voxel) / gradient)/r);
                        } else {
                            alpha = 0.;
                        }
                    }

                    totalAlpha *= 1. - alpha;

                    if (totalAlpha < 0.05) {
                        break;
                    }
                }

                voxelColor.r = tfEditor2D.triangleWidget.color.r;
                voxelColor.g = tfEditor2D.triangleWidget.color.g;
                voxelColor.b = tfEditor2D.triangleWidget.color.b;
                voxelColor.a = 1. - totalAlpha;

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                setPixels(i, j, resolutionScaling, depthScaling, pixelColor);
            }
        }
    }

    void transfer2DShade(double[] viewMatrix, int resolutionScaling, int depthScaling) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, (double)volume.getDimX() / 2., (double)volume.getDimY() / 2., (double)volume.getDimZ() / 2.);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        int viewDim = (int) Math.sqrt(Math.pow((double) volume.getDimZ(), 2.) +
                Math.pow((double) volume.getDimX(), 2.) +
                Math.pow((double) volume.getDimY(), 2.) );

        short voxel;
        VoxelGradient vGradient;
        double alpha, totalAlpha;
        TFColor current, composite;

        double r = tfEditor2D.triangleWidget.radius;
        double av = tfEditor2D.triangleWidget.color.a;
        short fv = tfEditor2D.triangleWidget.baseIntensity;

        for (int j = 0; j < image.getHeight(); j+=resolutionScaling) {
            for (int i = 0; i < image.getWidth(); i+=resolutionScaling) {

                totalAlpha = 1.;
                composite = new TFColor(0., 0., 0., 0.);
                //voxelColor.a = 0;

                for (int k = viewDim - 1; k >= 0; k-=depthScaling) {

                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];

                    voxel = getVoxelTriLinear(pixelCoord);
                    vGradient = getGradientTrilinear(pixelCoord);

                    if (vGradient.mag < tfEditor2D.triangleWidget.minGradMag || vGradient.mag > tfEditor2D.triangleWidget.maxGradMag) {
                        alpha = 0;
                    } else {
                        if (vGradient.mag == 0. && fv == voxel) {
                            alpha = av;
                        } else if (vGradient.mag > 0 && fv >= voxel - r * vGradient.mag && fv <= voxel + r * vGradient.mag) {
                            alpha = av * (1. - Math.abs((fv - voxel) / vGradient.mag)/r);
                        } else {
                            alpha = 0.;
                        }
                    }

                    current = shadePhong(tfEditor2D.triangleWidget.color, viewVec, vGradient, 0.1, 0.7, 0.2, 10.0);

                    composite.r = current.r * alpha + (1. - alpha) * composite.r;
                    composite.g = current.g * alpha + (1. - alpha) * composite.g;
                    composite.b = current.b * alpha + (1. - alpha) * composite.b;

                    totalAlpha *= 1. - alpha;
                }

                voxelColor.r = composite.r;
                voxelColor.g = composite.g;
                voxelColor.b = composite.b;
                voxelColor.a = 1. - totalAlpha;

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                setPixels(i, j, resolutionScaling, depthScaling, pixelColor);
            }
        }
    }

    void transfer2DShadeFront(double[] viewMatrix, int resolutionScaling, int depthScaling) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, (double)volume.getDimX() / 2., (double)volume.getDimY() / 2., (double)volume.getDimZ() / 2.);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        int viewDim = (int) Math.sqrt(Math.pow((double) volume.getDimZ(), 2.) +
                Math.pow((double) volume.getDimX(), 2.) +
                Math.pow((double) volume.getDimY(), 2.) );

        short voxel;
        VoxelGradient vGradient;
        double alpha;
        TFColor current, composite;

        double r = tfEditor2D.triangleWidget.radius;
        double av = tfEditor2D.triangleWidget.color.a;
        short fv = tfEditor2D.triangleWidget.baseIntensity;

        for (int j = 0; j < image.getHeight(); j+=resolutionScaling) {
            for (int i = 0; i < image.getWidth(); i+=resolutionScaling) {
                composite = new TFColor(0., 0., 0., 1.);

                for (int k = 0; k < viewDim; k += depthScaling) {

                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];

                    voxel = getVoxelTriLinear(pixelCoord);
                    vGradient = getGradientTrilinear(pixelCoord);

                    if (vGradient.mag < tfEditor2D.triangleWidget.minGradMag || vGradient.mag > tfEditor2D.triangleWidget.maxGradMag) {
                        alpha = 0;
                    } else {
                        if (vGradient.mag == 0. && fv == voxel) {
                            alpha = av;
                        } else if (vGradient.mag > 0 && fv >= voxel - r * vGradient.mag && fv <= voxel + r * vGradient.mag) {
                            alpha = av * (1. - Math.abs((fv - voxel) / vGradient.mag)/r);
                        } else {
                            alpha = 0.;
                        }
                    }

                    current = shadePhong(tfEditor2D.triangleWidget.color, viewVec, vGradient, 0.1, 0.7, 0.2, 10.0);

                    composite.r = composite.r + composite.a * alpha * current.r;
                    composite.g = composite.g + composite.a * alpha * current.g;
                    composite.b = composite.b + composite.a * alpha * current.b;
                    composite.a *= 1. - alpha;

                    if (composite.a < 0.05) {
                        break;
                    }
                }

                voxelColor.r = composite.r;
                voxelColor.g = composite.g;
                voxelColor.b = composite.b;
                voxelColor.a = 1. - composite.a;

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                setPixels(i, j, resolutionScaling, depthScaling, pixelColor);
            }
        }
    }

    TFColor shadePhong(TFColor color, double[] viewVec, VoxelGradient gradient, double ambient, double diff, double spec, double alpha) {
        double[] normal = new double[3];

        // normalization is not necessary since the view vector is already normalized
        //VectorMath.setVector(viewVec, -viewVec[0], -viewVec[1], -viewVec[2]);

        // In het simplified case that L=V we also have H=V, since V is normalized

        if (gradient.mag == 0.) {
            return new TFColor(ambient, ambient, ambient, 0.);
        } else {
            VectorMath.setVector(normal, gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag);
        }

        //double second = diff * VectorMath.dotproduct(viewVec, normal);
        double dotprod = VectorMath.dotproduct(normal, viewVec);

        if (/*second < 0 || */dotprod < 0) {
            return new TFColor(0., 0., 0., 0.);
        }

        double second = diff * dotprod;
        double third = spec * Math.pow(dotprod, alpha);

        return new TFColor(ambient + color.r * second + third, ambient + color.g * second + third, ambient + color.b * second + third, 0.0);
    }

    TFColor complexShadePhong(TFColor color, double[] viewVec, VoxelGradient gradient, double ambient, double diff, double spec, double alpha) {
        //double[] halfway = new double[3];
        double[] normal = new double[3];
        double[] reflect = new double[3];

        // normalization is not necessary since the view vector is already normalized
        //double viewVecLength = VectorMath.length(viewVec);

        //VectorMath.setVector(viewVec, -viewVec[0], -viewVec[1], -viewVec[2]);

        //VectorMath.setVector(halfway, viewVec[0] / viewVecLength, viewVec[1] / viewVecLength, viewVec[2] / viewVecLength);


        // In het simplified case that L=V we also have H=V, since V is normalized
        //VectorMath.setVector(halfway, viewVec[0], viewVec[1], viewVec[2]);

        if (gradient.mag == 0.) {

            VectorMath.setVector(normal, 0., 0., 0.);
            //return new TFColor(ambient, ambient, ambient, 0.);
        } else {
            VectorMath.setVector(normal, gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag);
        }

        //double second = diff * VectorMath.dotproduct(viewVec, normal);
        double dotprod = VectorMath.dotproduct(normal, viewVec);

        VectorMath.setVector(reflect, (2 * dotprod) * normal[0] - viewVec[0],
                (2 * dotprod) * normal[1] - viewVec[1],
                (2 * dotprod) * normal[2] - viewVec[2]);

        double dotprod2 = VectorMath.dotproduct(viewVec, reflect);

        if (dotprod2 < 0 || dotprod < 0) {
            return new TFColor(0., 0., 0., 0.);
        }

        double second = diff * dotprod;
        double third = spec * Math.pow(dotprod2, alpha);

        /*color.r = ambient + color.r * second + third;
        color.g = ambient + color.g * second + third;
        color.b = ambient + color.b * second + third;*/

        return new TFColor(ambient + color.r * ( second + third ),
                ambient + color.g * ( second + third ),
                ambient + color.b * ( second + third ), 0.0);
    }

    private void setPixels(int i, int j, int resolutionScaling, int depthScaling, int pixelColor) {
        if (resolutionScaling > 1) {
            int odd = resolutionScaling % 2;

            if (odd == 0) {
                int half = resolutionScaling/2;

                for (int i_k = -half + 1; i_k <= half; ++i_k) {
                    for (int j_k = -half + 1; j_k <= half; ++j_k) {
                        if ((i + i_k >= 0) && (i + i_k < image.getWidth()) &&
                                (j + j_k >= 0) && (j + j_k) < image.getHeight()) {

                            image.setRGB(i + i_k, j + j_k, pixelColor);
                        }
                    }
                }
            } else {
                int half = (resolutionScaling-1)/2;

                for (int i_k = -half; i_k <= half; ++i_k) {
                    for (int j_k = -half; j_k <= half; ++j_k) {
                        if ((i + i_k >= 0) && (i + i_k < image.getWidth()) &&
                                (j + j_k >= 0) && (j + j_k) < image.getHeight()) {

                            image.setRGB(i + i_k, j + j_k, pixelColor);
                        }
                    }
                }
            }
        } else {
            image.setRGB(i, j, pixelColor);
        }
    }

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();

        int resolutionScaling = 1, depthScaling = 1;

        if (this.interactiveMode) {
            resolutionScaling = 3;
            depthScaling = 4;
        }

        if (panel.mode == 0) {

            if (panel.shade) {
                colorSlicer(viewMatrix);
            } else {
                slicer(viewMatrix);
            }
        } else if (panel.mode == 1) {
            mip(viewMatrix, resolutionScaling, depthScaling);
        } else if (panel.mode == 2) {
            //composite(viewMatrix, resolutionScaling, depthScaling);
            compositeFront(viewMatrix, resolutionScaling, depthScaling);
        } else if (panel.mode == 3) {

            if (panel.shade) {
                //transfer2DShade(viewMatrix, resolutionScaling, depthScaling);
                transfer2DShadeFront(viewMatrix, resolutionScaling, depthScaling);
            } else {
                transfer2D(viewMatrix, resolutionScaling, depthScaling);
            }
        }

        // uncomment to save a png of the generated image (out.png)
        /*
        try {
            ImageIO.write(image, "png",new File("out.png"));
        } catch (java.io.IOException e) {
            System.out.println("Could not save image to file.");
            System.out.println(e.getMessage());
        }
        */


        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
