// Triangle class
// defines a Triangle shape

import javax.vecmath.*;

public class Triangle extends Shape {
	private Vector3f p0, p1, p2;	// three vertices make a triangle
	private Vector3f n0, n1, n2;	// normal at each vertex

	public Triangle() {
	}
	public Triangle(Vector3f _p0, Vector3f _p1, Vector3f _p2, Material mat) {
		p0 = new Vector3f(_p0);
		p1 = new Vector3f(_p1);
		p2 = new Vector3f(_p2);
		material = mat;
		Vector3f normal = new Vector3f();
		Vector3f v1 = new Vector3f();
		Vector3f v2 = new Vector3f();
		v1.sub(p1, p0);
		v2.sub(p2, p0);
		normal.cross(v1, v2);
		normal.normalize();				// compute default normal:
		n0 = new Vector3f(normal);		// the normal of the plane defined by the triangle
		n1 = new Vector3f(normal);
		n2 = new Vector3f(normal);
	}
	public Triangle(Vector3f _p0, Vector3f _p1, Vector3f _p2,
					Vector3f _n0, Vector3f _n1, Vector3f _n2,
					Material mat) {
		p0 = new Vector3f(_p0);
		p1 = new Vector3f(_p1);
		p2 = new Vector3f(_p2);
		material = mat;
		n0 = new Vector3f(_n0);		// the normal of the plane defined by the triangle
		n1 = new Vector3f(_n1);
		n2 = new Vector3f(_n2);
	}
	public HitRecord hit(Ray ray, float tmin, float tmax) {

		/* YOUR WORK HERE: complete the triangle's intersection routine
		 * Normal should be computed by a bilinear interpolation from n0, n1 and n2
		 * using the barycentric coordinates: alpha, beta, (1.0 - alpha - beta) */
		float alpha, beta, gamma;
		
//		Vector3f p = new Vector3f((p0.x+p1.x+p2.x)/3, (p0.y+p1.y+p2.y)/3, (p0.z+p1.z+p2.z)/3);
//		Vector3f v0 = new Vector3f(p1);
//		p1.sub(p2);
//		Vector3f v1 = new Vector3f(p2);
//		p2.sub(p0);
//		Vector3f v2 = new Vector3f(p);
//		p.sub(p0);
//	    float d00 = v0.dot(v0);
//	    float d01 = v0.dot(v1);
//	    float d11 = v1.dot(v1);
//	    float d20 = v2.dot(v0);
//	    float d21 = v2.dot(v1);
//	    float denom = d00 * d11 - d01 * d01;
//	    alpha = (d11 * d20 - d01 * d21) / denom;
//	    beta = (d00 * d21 - d01 * d20) / denom;
//	    gamma = 1.0f - alpha - beta;
		
		HitRecord hitRec = null;

		//COLUMN 1
		Vector3f dir = ray.getDirection();
		//COLUMN 2
		Vector3f p2minusp0 = new Vector3f(p2);
		p2minusp0.sub(p0);
		//COLUMN 3
		Vector3f p2minusp1 = new Vector3f(p2);
		p2minusp1.sub(p1);
		//COLUMN 4
		Vector3f p2minusOrigin = new Vector3f(p2);
		p2minusOrigin.sub(ray.getOrigin());
		
		float[][] matrixDen = new float[3][3];
		matrixDen[0][0]=dir.x; matrixDen[0][1]=p2minusp0.x; matrixDen[0][2]=p2minusp1.x;
		matrixDen[1][0]=dir.y; matrixDen[1][1]=p2minusp0.y; matrixDen[1][2]=p2minusp1.y;
		matrixDen[2][0]=dir.z; matrixDen[2][1]=p2minusp0.z; matrixDen[2][2]=p2minusp1.z;
		float[][] matrix_xNum = new float[3][3];
		matrix_xNum[0][0]=p2minusOrigin.x; matrix_xNum[0][1]=p2minusp0.x; matrix_xNum[0][2]=p2minusp1.x;
		matrix_xNum[1][0]=p2minusOrigin.y; matrix_xNum[1][1]=p2minusp0.y; matrix_xNum[1][2]=p2minusp1.y;
		matrix_xNum[2][0]=p2minusOrigin.z; matrix_xNum[2][1]=p2minusp0.z; matrix_xNum[2][2]=p2minusp1.z;
		float[][] matrix_yNum = new float[3][3];
		matrix_yNum[0][0]=dir.x; matrix_yNum[0][1]=p2minusOrigin.x; matrix_yNum[0][2]=p2minusp1.x;
		matrix_yNum[1][0]=dir.y; matrix_yNum[1][1]=p2minusOrigin.y; matrix_yNum[1][2]=p2minusp1.y;
		matrix_yNum[2][0]=dir.z; matrix_yNum[2][1]=p2minusOrigin.z; matrix_yNum[2][2]=p2minusp1.z;
		float[][] matrix_zNum = new float[3][3];
		matrix_zNum[0][0]=dir.x; matrix_zNum[0][1]=p2minusp0.x; matrix_zNum[0][2]=p2minusOrigin.x;
		matrix_zNum[1][0]=dir.y; matrix_zNum[1][1]=p2minusp0.y; matrix_zNum[1][2]=p2minusOrigin.y;
		matrix_zNum[2][0]=dir.z; matrix_zNum[2][1]=p2minusp0.z; matrix_zNum[2][2]=p2minusOrigin.z;
		
		float det_MatrixDen = determinant(matrixDen);
		float det_matrix_xNum = determinant(matrix_xNum);
		float det_matrix_yNum = determinant(matrix_yNum);
		float det_matrix_zNum = determinant(matrix_zNum);
		
		float t_x = det_matrix_xNum/det_MatrixDen;
		float alpha_y = det_matrix_yNum/det_MatrixDen;
		float beta_z = det_matrix_zNum/det_MatrixDen;
		
		alpha = alpha_y;
		beta = beta_z;
		gamma = (1 - alpha - beta);
		
		hitRec = new HitRecord();
		//Determinant != 0 & ABG >= 0
		if(det_MatrixDen != 0 && alpha >=0 && beta >=0 && t_x >=0 && (alpha+beta) < 1){
			Vector3f intPoint = ray.pointAt(t_x);
			//T
			hitRec.t = t_x;
			//POSITION
			hitRec.pos = intPoint;
			//NORMAL
			Vector3f pNormal = new Vector3f((alpha*n0.x + beta*n1.x + gamma*n2.x),
											(alpha*n0.y + beta*n1.y + gamma*n2.y),
											(alpha*n0.z + beta*n1.z + gamma*n2.z));    	
	    	
	    	hitRec.normal = pNormal;
	    	//MATERIAL
	    	hitRec.material = material;
		}
		
		return hitRec;
	}
	public static float determinant(float[][] matrix){ //method sig. takes a matrix (two dimensional array), returns determinant.
	    float sum=0; 
	    int s;
	    if(matrix.length==1){  //bottom case of recursion. size 1 matrix determinant is itself.
	      return(matrix[0][0]);
	    }
	    for(int i=0;i<matrix.length;i++){ //finds determinant using row-by-row expansion
	      float[][]smaller= new float[matrix.length-1][matrix.length-1]; //creates smaller matrix- values not in same row, column
	      for(int a=1;a<matrix.length;a++){
	        for(int b=0;b<matrix.length;b++){
	          if(b<i){
	            smaller[a-1][b]=matrix[a][b];
	          }
	          else if(b>i){
	            smaller[a-1][b-1]=matrix[a][b];
	          }
	        }
	      }
	      if(i%2==0){ //sign changes based on i
	        s=1;
	      }
	      else{
	        s=-1;
	      }
	      sum+=s*matrix[0][i]*(determinant(smaller)); //recursive step: determinant of larger determined by smaller.
	    }
	    return(sum); //returns determinant value. once stack is finished, returns final determinant.
	  }
	//void Barycentric(Point p, Point a, Point b, Point c, float &u, float &v, float &w)
}
