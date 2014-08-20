varying vec3 normal;
varying vec4 vertex;
varying vec3 halfvec;

void prepareLight()
{
	/* Calculate the normal */
	normal = normalize(gl_NormalMatrix * gl_Normal).xyz;

	/* Transform the vertex position to eye space */
	vertex = vec4(gl_ModelViewMatrix * gl_Vertex);

	vec3 L = normalize(gl_LightSource[0].position.xyz - vertex.xyz);

	halfvec = vec3(L + vec3(0.0, 0.0, 1.0));
}

