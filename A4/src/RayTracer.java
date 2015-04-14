// RayTracer class

import javax.vecmath.*;

import java.util.*;
import java.awt.*;
import java.awt.image.*;

import javax.imageio.*;

import java.io.*;

public class RayTracer {

	private Color3f image[][];	// image that stores floating point color
	private String image_name;	// output image name
	private int width, height;	// image width, height
	private int xsample, ysample;	// samples used for super sampling
	private Color3f background;	// background color
	private Color3f ambient;	// ambient color
	private int maxdepth;		// max recursion depth for recursive ray tracing
	private float exposure;		// camera exposure for the entire scene

	private Camera camera;
	private Vector<Material> materials = new Vector<Material> ();	// array of materials
	private Vector<Shape> shapes = new Vector<Shape> ();			// array of shapes
	private Vector<Light> lights = new Vector<Light> ();			// array of lights

	private void initialize() {
		width = 256;
		height = 256;
		xsample = 1;
		ysample = 1;
		maxdepth = 5;
		background = new Color3f(0,0,0);
		ambient = new Color3f(0,0,0);
		exposure = 1.0f;

		image_name = new String("output.png");

		camera = new Camera(new Vector3f(0,0,0), new Vector3f(0,-1,0), new Vector3f(0,1,0), 45.f, 1.f);

		// add a default material: diffuse material with constant 1 reflectance
		materials.add(Material.makeDiffuse(new Color3f(0,0,0), new Color3f(1,1,1)));
	}

	public static void main(String[] args) {
		if (args.length == 1) {
			new RayTracer(args[0]);
		} else {
			System.out.println("Usage: java RayTracer input.scene");
		}
	}

	/*CALLED FOR EACH RAY*/
	private Color3f raytracing(Ray ray, int depth)
	{
		/* YOUR WORK HERE: complete the ray tracing function
		 * Feel free to make new function as needed. For example, you may want to add a 'shading' function */
		Color3f returnColor = new Color3f(0,0,0);

		//SET HIT = WHERE RAY INTERSECTS SHAPE
		HitRecord hit = IntersectScene(ray, shapes);
//		if(depth == 1){
//			if(hit==null)
//				System.out.println("Hit is NULL");
//			else
//				System.out.println(hit.t);
//		}
		if(hit != null){
			//REFLECT
			Vector3f light_pos = new Vector3f();
			Vector3f light_dir = new Vector3f();
			lights.elementAt(0).getLight(hit.pos, light_pos, light_dir);

			Ray sr = new Ray(hit.pos, light_dir);
			Vector3f tempLightPos = new Vector3f(light_pos);
			tempLightPos.sub(hit.pos);
			float light_dist = tempLightPos.length();
			
			//System.out.println(hit.material);
			Color3f materialKr = new Color3f(hit.material.Kr);
			Color3f materialKt = new Color3f(hit.material.Kt);
			Color3f returnReflRefr = new Color3f(0,0,0);

			Vector3f L = new Vector3f(ray.d);
			L.scale(-1);
			Vector3f eye = new Vector3f(ray.d);
			Vector3f N = new Vector3f(hit.normal);
			Vector3f V = new Vector3f(ray.getDirection());
			
			Vector3f tempL = new Vector3f(light_dir);
			Vector3f tempN = new Vector3f(hit.normal);
			float tempR = tempL.dot(tempN);
			tempR*=2;
			tempN.scale(tempR);
			tempN.sub(tempL);
			Vector3f R = new Vector3f(tempN);
			
			Vector3f nL = new Vector3f(eye);
			
			

			int maxDepth = 5;
			if ((materialKr.x+materialKr.y+materialKr.z >= 0.0001f) && depth < maxDepth){
				Vector3f refl = reflect(L, N);
				//System.out.println("nL="+nL+"  N="+N+"  ior="+hit.material.ior);
				Ray refl_ray = new Ray(hit.pos, refl);
				//System.out.println(hit.material.ior);
				//System.out.println(hit.material.ior);
				//System.out.println(refr);				
				returnReflRefr.x += materialKr.x * raytracing(refl_ray, depth+1).x;
				
				returnReflRefr.y += materialKr.y * raytracing(refl_ray, depth+1).y;
				returnReflRefr.z += materialKr.z * raytracing(refl_ray, depth+1).z;
				returnColor = new Color3f(returnReflRefr.x, returnReflRefr.y, returnReflRefr.z);
			}
			else if ((materialKt.x+materialKt.y+materialKt.z >= 0.0001f) && depth < maxDepth){
				Vector3f refr = refract(nL, N, hit.material.ior);
				//System.out.println(hit.material.ior);
				//System.out.println(hit.material.ior);
				//System.out.println(refr);
				Ray refr_ray = new Ray(hit.pos, refr);
				returnReflRefr.x += materialKt.x * raytracing(refr_ray, depth+1).x;
				if(depth == 0) {
					depth += 0;
				}
				returnReflRefr.y += materialKt.y * raytracing(refr_ray, depth+1).y;
				returnReflRefr.z += materialKt.z * raytracing(refr_ray, depth+1).z;
				returnColor = new Color3f(returnReflRefr.x, returnReflRefr.y, returnReflRefr.z);
			}
			else{
				//RETURN COLOR OF WHAT RAY HIT
				returnColor = getShadingColor(ray, hit);
			}
		}
		else
			returnColor = background;
		/*if(hit!=null)
			return new Color3f(hit.normal.x, hit.normal.y, hit.normal.z);
		else return new Color3f(0,0,0);*/
		return returnColor;
	}
	
	private HitRecord IntersectScene(Ray ray, Vector<Shape> shapesInScene){
		
		HitRecord hit = new HitRecord();
		float tmax = Float.MAX_VALUE;
		int number_of_shapes = shapesInScene.size();
		hit = null;
		for (int i=0; i < number_of_shapes; i++){
			HitRecord shape_hit = null;
			//HitRecord temp_hit;
								//SPHERE.HIT
			shape_hit = shapesInScene.elementAt(i).hit(ray, 0.0001f, tmax);
			//temp_hit = shapesInScene.elementAt(i);
			float tmin = 0.00001f;
					if (shape_hit != null && shape_hit.t < tmax && shape_hit.t > tmin){
						tmax = shape_hit.t;
						hit = new HitRecord();
						hit.set(shape_hit);
						//System.out.println("HIT!");
					}
		}
		return hit;
	}
	//If intersect returnColor with shading & specular
	private Color3f getShadingColor(Ray ray, HitRecord hit){
		Color3f c = new Color3f(0,0,0);
		
		//AMBIENT REFLECTANCE
		Color3f materialKa = new Color3f(hit.material.Ka);

		//I * ka
		float materialKaRed = materialKa.x;
		float materialKaGreen = materialKa.y;
		float materialKaBlue = materialKa.z;
		Color3f B_Refl_Ka = new Color3f(ambient.x * materialKaRed,
										ambient.y * materialKaGreen,
										ambient.z * materialKaBlue);
		
		//FOR EACH LIGHT
		for (int i = 0; i < lights.size(); i ++) {
			Vector3f light_pos = new Vector3f();
			Vector3f light_dir = new Vector3f();
			//System.out.println(lights.elementAt(i).getLight(hit.pos, light_pos, light_dir));
			lights.elementAt(i).getLight(hit.pos, light_pos, light_dir);
			
			Ray sr = new Ray(hit.pos, light_dir);
			Vector3f tempLightPos = new Vector3f(light_pos);
			tempLightPos.sub(hit.pos);
			float light_dist = tempLightPos.length();
			HitRecord shadow_hit = IntersectScene(sr, shapes);
			
			if (shadow_hit != null && shadow_hit.t <= light_dist)
				continue;
			else{
				 Vector3f L = new Vector3f(light_dir);
				 Vector3f N = new Vector3f(hit.normal);
				 Vector3f V = new Vector3f(ray.getDirection());
				 V.scale(-1);
				 Vector3f tempL = new Vector3f(light_dir);
				 Vector3f tempN = new Vector3f(hit.normal);
				 float tempR = tempL.dot(tempN);
				 tempR*=2;
				 tempN.scale(tempR);
				 tempN.sub(tempL);
				 Vector3f R = new Vector3f(tempN);
				 
				 
				 
				 float lightRed = lights.elementAt(i).intensity.x / 255;
				 float lightGreen = lights.elementAt(i).intensity.y / 255;
				 float lightBlue = lights.elementAt(i).intensity.z / 255;
				 				 
				 //DIFFUSE REFLECTANCE
				 Color3f materialKd = new Color3f(hit.material.Kd);
				 //System.out.println(hit.material.Kd);
				 //I * ka * max(dot(N,L),0) 
				 float materialKdRed = materialKd.x; //Materials Diffuse Reflectance
				 float materialKdGreen = materialKd.y;
				 float materialKdBlue = materialKd.z;
				 float test5 = Math.max(N.dot(L), 0);
				 float B_Refl_KdR = lightRed * materialKdRed * Math.max(N.dot(L), 0);
				 float B_Refl_KdG = lightGreen * materialKdGreen * Math.max(N.dot(L), 0);
				 float B_Refl_KdB = lightBlue * materialKdBlue * Math.max(N.dot(L), 0); 
				 Color3f B_Refl_Kd = new Color3f(B_Refl_KdR, B_Refl_KdG, B_Refl_KdB);
				 
				 //SPECULAR REFLECTANCE
				 Color3f materialKs = new Color3f(hit.material.Ks);
				 //I * ka * max(dot(N,L),0) 
				 float materialKsRed = materialKs.x;
				 float materialKsGreen = materialKs.y;
				 float materialKsBlue = materialKs.z;
				 Color3f iXks = new Color3f(lightRed * materialKsRed,
						 					lightGreen * materialKsGreen,
						 					lightBlue * materialKsBlue);
				 float test2 = R.dot(V);
				 float test = (float)Math.max(R.dot(V), 0);
				 float B_Spec_Ks = (float)Math.pow(Math.max(R.dot(V), 0), hit.material.phong_exp);
				 Color3f B_Refl_Ks = new Color3f(iXks.x * B_Spec_Ks,
						 					iXks.y * B_Spec_Ks,
						 					iXks.z * B_Spec_Ks);
				 //System.out.println(hit.material.phong_exp);
				 
				  c = new Color3f((B_Refl_Kd.x + B_Refl_Ka.x + B_Refl_Ks.x), (B_Refl_Kd.y + B_Refl_Ka.y + B_Refl_Ks.y), (B_Refl_Kd.z + B_Refl_Ka.z + B_Refl_Ks.z));
				 //c = new Color3f(hit.material.Kd);
				 //System.out.println(hit.material.Kd.x);
				 //c = new Color3f(B_Refl_Kd.x, B_Refl_Kd.y, B_Refl_Kd.z);
				 //c = new Color3f(B_Refl_Ka.x, B_Refl_Ka.y, B_Refl_Ka.z);
				 //c = new Color3f(B_Refl_Ks.x, B_Refl_Ks.y, B_Refl_Ks.z);
				 //c = new Color3f((B_Refl_Kd.x + B_Refl_Ka.x), (B_Refl_Kd.y + B_Refl_Ka.y), (B_Refl_Kd.z + B_Refl_Ka.z));
			}
		}
		//System.out.println(c);
		return c;
	}
	

	// reflect a direction (in) around a given normal
	/* NOTE: dir is assuming to point AWAY from the hit point
	 * if your ray direction is point INTO the hit point, you should flip
	 * the sign of the direction before calling reflect
	 */
	private Vector3f reflect(Vector3f dir, Vector3f normal)
	{
		Vector3f out = new Vector3f(normal);
		out.scale(2.f * dir.dot(normal));
		out.sub(dir);
		return out;
	}

	// refract a direction (in) around a given normal and 'index of refraction' (ior)
	/* NOTE: dir is assuming to point INTO the hit point
	 * (this is different from the reflect function above, which assumes dir is pointing away
	 */
	private Vector3f refract(Vector3f dir, Vector3f normal, float ior)
	{
		float mu;
		mu = (normal.dot(dir) < 0) ? 1.f / ior : ior;

		float cos_thetai = dir.dot(normal);
		float sin_thetai2 = 1.f - cos_thetai*cos_thetai;

		//System.out.println("st = "+sin_thetai2);
		
		float sin_thetar = mu*(float)Math.sqrt(sin_thetai2);
		float cos_thetar = (float)Math.sqrt(1.f - sin_thetar*sin_thetar);

		if (mu*mu*sin_thetai2 > 1.f){
			//System.out.println("----------- = "+(mu*mu*sin_thetai2) + ", " + (normal.dot(dir)) + ", " + mu);
			//mu = 1.f / mu;
			return null;
		}
		Vector3f out = new Vector3f(normal);
		if (cos_thetai > 0)
		{
			out.scale(-mu*cos_thetai+cos_thetar);
			out.scaleAdd(mu, dir, out);

		} else {

			out.scale(-mu*cos_thetai-cos_thetar);
			out.scaleAdd(mu, dir, out);
		}
		out.normalize();
		return out;
	}

	public RayTracer(String scene_name) {

		// initialize and set default parameters
		initialize();

		// parse scene file
		parseScene(scene_name);

		// create floating point image
		image = new Color3f[width][height];

		int i, j;
		float x, y;
		for (j=0; j<height; j++)
		//for(j = 162; j<=162; j++)	
		{
			y = (float)j / (float)height;
			System.out.print("\rray tracing... " + j*100/height + "%");
			
			for (i=0; i<width; i ++)
			//for(i = 248; i<=248; i++)
			{
				// System.out.println("uv: " +i + ", " + j);
				x = (float)i / (float)width;
				image[i][j] = raytracing(camera.getCameraRay(x, y), 0);
				int z=0;
			}
		}
		System.out.println("\rray tracing completed.                       ");
				
		writeImage();
	}

	private void parseScene(String scene_name)
	{
		File file = null;
		Scanner scanner = null;
		try {
			file = new File(scene_name);
			scanner = new Scanner(file);
		} catch (IOException e) {
			System.out.println("error reading from file " + scene_name);
			System.exit(0);
		}
		String keyword;
		while(scanner.hasNext()) {

			keyword = scanner.next();
			// skip the comment lines
			if (keyword.charAt(0) == '#') {
				scanner.nextLine();
				continue;
			}
			if (keyword.compareToIgnoreCase("image")==0) {

				image_name = scanner.next();
				width = scanner.nextInt();
				height = scanner.nextInt();
				exposure = scanner.nextFloat();

			} else if (keyword.compareToIgnoreCase("camera")==0) {

				Vector3f eye = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
				Vector3f at  = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
				Vector3f up  = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
				float fovy = scanner.nextFloat();
				float aspect_ratio = (float)width / (float)height;

				camera = new Camera(eye, at, up, fovy, aspect_ratio);

			} else if (keyword.compareToIgnoreCase("background")==0) {

				background = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());

			} else if (keyword.compareToIgnoreCase("ambient")==0) { 

				ambient = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());

			} else if (keyword.compareToIgnoreCase("maxdepth")==0) {

				maxdepth = scanner.nextInt();

			} else if (keyword.compareToIgnoreCase("light")==0) {

				// parse light
				parseLight(scanner);

			} else if (keyword.compareToIgnoreCase("material")==0) {

				// parse material
				parseMaterial(scanner);

			} else if (keyword.compareToIgnoreCase("shape")==0) {

				// parse shape
				parseShape(scanner);
		
			} else {
				System.out.println("undefined keyword: " + keyword);
			}
		}
		scanner.close();
	}

	private void parseLight(Scanner scanner)
	{
		String lighttype;
		lighttype = scanner.next();
		if (lighttype.compareToIgnoreCase("point")==0) {

			/* add a new point light */
			Vector3f pos = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Color3f intens = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			lights.add(new PointLight(pos, intens));

		} else if (lighttype.compareToIgnoreCase("spot")==0) {

			/* add a new spot light */
			Vector3f from = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Vector3f to = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			float spot_exponent = scanner.nextFloat();
			float spot_cutoff = scanner.nextFloat();
			Color3f intens = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());

			lights.add(new SpotLight(from, to, spot_exponent, spot_cutoff, intens));

		} else if (lighttype.compareToIgnoreCase("area")==0) {

			/* YOUR WORK HERE: complete the area light
			 * Note that you do not need to create a new type of light source.
			 * Instead, you will convert an area light
			 * to a collection of point lights and add them all to the 'lights' array 
			 * 
			 * ADd all 256 point lights here
			 * Compute pos & intensit
			 * 
			 * */

		} else {
			System.out.println("undefined light type: " + lighttype);
		}
	}

	private void parseMaterial(Scanner scanner)
	{
		String mattype;
		mattype = scanner.next();
		if (mattype.compareToIgnoreCase("diffuse")==0) {

			Color3f ka = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Color3f kd = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			materials.add(Material.makeDiffuse(ka, kd));

		} else if (mattype.compareToIgnoreCase("specular")==0) {

			Color3f ka = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Color3f kd = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Color3f ks = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			float phong_exp = scanner.nextFloat();
			materials.add(Material.makeSpecular(ka, kd, ks, phong_exp));

		} else if (mattype.compareToIgnoreCase("mirror")==0) {

			Color3f kr = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			materials.add(Material.makeMirror(kr));

		} else if (mattype.compareToIgnoreCase("glass")==0) {

			Color3f kr = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Color3f kt = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			float ior = scanner.nextFloat();
			materials.add(Material.makeGlass(kr, kt, ior));

		} else if (mattype.compareToIgnoreCase("super")==0) {

			Color3f ka = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Color3f kd = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Color3f ks = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			float phong_exp = scanner.nextFloat();
			Color3f kr = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Color3f kt = new Color3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			float ior = scanner.nextFloat();
			materials.add(Material.makeSuper(ka, kd, ks, phong_exp, kr, kt, ior));			
		}

		else {
			System.out.println("undefined material type: " + mattype);
		}

	}

	private void parseShape(Scanner scanner)
	{
		String shapetype;
		shapetype = scanner.next();
		Material material = materials.lastElement();
		if (shapetype.compareToIgnoreCase("plane")==0) {

			Vector3f P0 = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Vector3f N = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			shapes.add(new Plane(P0, N, material));

		} else if (shapetype.compareToIgnoreCase("sphere")==0) {

			Vector3f center = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			float radius = scanner.nextFloat();
			shapes.add(new Sphere(center, radius, material));

		} else if (shapetype.compareToIgnoreCase("triangle")==0) {

			Vector3f p0 = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Vector3f p1 = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Vector3f p2 = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			shapes.add(new Triangle(p0, p1, p2, material));

		} else if (shapetype.compareToIgnoreCase("triangle_n")==0) {

			Vector3f p0 = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Vector3f p1 = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Vector3f p2 = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());

			Vector3f n0 = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Vector3f n1 = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());
			Vector3f n2 = new Vector3f(scanner.nextFloat(), scanner.nextFloat(), scanner.nextFloat());

			shapes.add(new Triangle(p0, p1, p2, n0, n1, n2, material));

		} else if (shapetype.compareToIgnoreCase("trimesh")==0) {

			TriMesh	mesh = new TriMesh();
			mesh.load(scanner.next());

			if (mesh.type.compareToIgnoreCase("triangle")==0) {
				int i;
				int idx0, idx1, idx2;
				for (i=0; i<mesh.faces.length/3; i++) {
					idx0 = mesh.faces[i*3+0];
					idx1 = mesh.faces[i*3+1];
					idx2 = mesh.faces[i*3+2];
					shapes.add(new Triangle(mesh.verts[idx0], mesh.verts[idx1], mesh.verts[idx2], material));
				}

			} else if (mesh.type.compareToIgnoreCase("triangle_n")==0) {
				int i;
				int idx0, idx1, idx2;
				for (i=0; i<mesh.faces.length/3; i++) {
					idx0 = mesh.faces[i*3+0];
					idx1 = mesh.faces[i*3+1];
					idx2 = mesh.faces[i*3+2];
					shapes.add(new Triangle(mesh.verts[idx0], mesh.verts[idx1], mesh.verts[idx2],
											mesh.normals[idx0], mesh.normals[idx1], mesh.normals[idx2],
											material));
				}

			} else {
				System.out.println("undefined trimesh type: " + mesh.type);
			}


		} else {
			System.out.println("undefined shape type: " + shapetype);
		}
	}

	// write image to a disk file
	// image will be multiplied by exposure
	private void writeImage() {
		int x, y, index;
		int pixels[] = new int[width * height];

		index = 0;
		// apply a standard 2.2 gamma correction
		float gamma = 1.f / 2.2f;
		for (y=height-1; y >= 0; y --) {
			for (x=0; x<width; x ++) {
				Color3f c = new Color3f(image[x][y]);
				c.x = (float)Math.pow(c.x*exposure, gamma);
				c.y = (float)Math.pow(c.y*exposure, gamma);
				c.z = (float)Math.pow(c.z*exposure, gamma);
				c.clampMax(1.f);
				pixels[index++] = c.get().getRGB();

			}
		}

		BufferedImage oimage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		oimage.setRGB(0, 0, width, height, pixels, 0, width);
		File outfile = new File(image_name);
		try {
			ImageIO.write(oimage, "png", outfile);
		} catch(IOException e) {
		}
	}
}
