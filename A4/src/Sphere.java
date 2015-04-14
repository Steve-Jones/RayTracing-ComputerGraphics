// Sphere class
// defines a Sphere shape

import java.util.ArrayList;

import javax.vecmath.*;

public class Sphere extends Shape {
	private Vector3f center;	// center of sphere
	private float radius;		// radius of sphere

	public Sphere() {
	}
	public Sphere(Vector3f pos, float r, Material mat) {
		center = new Vector3f(pos);
		radius = r;
		material = mat;
	}
	public HitRecord hit(Ray ray, float tmin, float tmax) {
		//ray.Origin - center) * ray.direction
		/* YOUR WORK HERE: complete the sphere's intersection routine */
		
		//Find A
		float A = (float)(ray.getDirection().lengthSquared());
		//System.out.println("DIRECTION="+ray.getDirection().x);
		//Find B
		Vector3f tempVectorB = new Vector3f(ray.getOrigin());
		tempVectorB.sub(center);
		tempVectorB.scale(2);
		float B = tempVectorB.dot(ray.getDirection());
		//Find C
		Vector3f tempVectorC = new Vector3f(ray.getOrigin());
		tempVectorC.sub(center);
		float C = (float)tempVectorC.lengthSquared() - (float)(radius*radius);
		//Quadratic Eq.
		float t1;
		float t2;
//		System.out.println("A="+A+"  B="+B+"  C="+C);
//		System.out.println("NAN = "+(B*B-(4*A*C)));
		if((B*B-(4*A*C)) < 0){
			t1 = -1;
			t2 = -1;
		}
		else{
			t1 = (float) ((-B + Math.sqrt(B*B-(4*A*C)))/2.0*A);
			t2 = (float) ((-B - Math.sqrt(B*B-(4*A*C)))/2.0*A);
		}
	    
	    HitRecord hitRec = null;
	    float returnT=10000;
	    float closeT = 1000000;
	    int i=0;
	    
	    //If t1 is closer and not negative
	    if(t1<closeT && tmin <= t1 && t1 <= tmax){
	    	closeT=t1;
	    }
	    if(t2<closeT && tmin <= t2 && t2 <= tmax){
	    	closeT=t2;
	    }
	    if(closeT < 999999){
		    returnT = closeT;
	    	hitRec = new HitRecord();
	    	//Set t
	    	Vector3f intPoint = ray.pointAt(returnT);
	    	hitRec.t = returnT;
	    	//Set pos
	    	hitRec.pos = ray.pointAt(returnT);
	    	//Set normal
	    	intPoint.sub(center);
	    	Vector3f returnTnormal = new Vector3f(intPoint);
	    	returnTnormal.normalize();
	    	hitRec.normal = returnTnormal;
	    	//Set material
	    	hitRec.material = material;
		    //System.out.println("\nt1 = "+t1+"\nt2 = "+t2);
		    //RETURN hitRec
	    }
	    return hitRec;
	}
}
