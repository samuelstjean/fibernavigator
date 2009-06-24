///////////////////////////////////////////////////////////////////////////////
 //
 ///// /////  ////       Project  :   FAnToM
 //   //  // //  //      Module   :   Fge (Rendering and Viewer Components)
 //   //  // /////       File     :   $RCSfile: $
 //   //  // //          Language :   C++
 //    /////  ////       Date     :   $Date: $
 //             Author   :   $Author: ebaum $
 //////////              Revision :   $Revision: 8037 $

varying vec3 TexCoord;

uniform int dimX, dimY, dimZ;

uniform sampler3D texes[10];

uniform float threshold[10];
uniform int type[10];
uniform int countTextures;

void lookupTex(inout vec4 color, in int type, in sampler3D tex, in float threshold, in vec3 v)
{
	vec3 col1;

	if (type == 1)
	{
		col1.r = clamp( texture3D(tex, v).r, 0.0, 1.0);
		col1.g = clamp( texture3D(tex, v).g, 0.0, 1.0);
		col1.b = clamp( texture3D(tex, v).b, 0.0, 1.0);

		if ( (length(col1) - threshold) < 0.0)
		{
			discard;
		}
	}

	if (type == 4)
	{

	}

}

/////////////////////////////////////////////////////////////////////////////////////////////
 // Transformation -- fragment shader -- main
 //
 // Gets called for every fragment and uses color and transforms the tennsor to image space
 /////////////////////////////////////////////////////////////////////////////////////////////
void main() {
	vec4 color = vec4(0.0);

	vec3 v = TexCoord;
	v.x = (v.x) / float( dimX);
	v.y = (v.y) / float( dimY);
	v.z = (v.z) / float( dimZ);

	for (int i = 9; i > -1; i--) {
		lookupTex(color, type[i], texes[i], threshold[i], v);
	}

	gl_FragColor = gl_Color;
}
