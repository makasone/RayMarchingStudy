Texture2D<float4>	tex0		: register(t0);
//Texture2D<float4>	tex1		: register(t1);
TextureCube texCube : register(t1);

SamplerState		sampler0	: register(s0);


//-------------------------------------------------------------------------------------
cbuffer System : register(b0)
{
	float g_time;

}

struct VSInput
{
	float4	pos : POSITION;
	float2	uv	: TEXCOORD0;
};

struct PSInput {
	float4	position	: SV_POSITION;
	float2	uv			: TEXCOORD0;
};

//----------------------------------------------------------------------------------------------
#define PI	3.14159265359
#define PI2	PI * 2.0
#define XAxis float3(1.0, 0.0, 0.0)
#define YAxis float3(0.0, 1.0, 0.0)
#define ZAxis float3(0.0, 0.0, 1.0)

//----------------------------------------------------------------------------------------------
struct Ray
{
	float3 start;
	float3 dir;
};

struct HitInfo
{
	Ray ray;
	float3 pos;
	float3 normal;
	float distance;

	float material;
};

struct DirectionalLight
{
	float3 direction;
	float3 color;
	float scale;
};

struct Scene
{
	DirectionalLight dirLight;
};

//----------------------------------------------------------------------------------------------
float3x3 calculateEyeRayTransformationMatrix(in float3 ro, in float3 ta, in float roll)
{
	float3 ww = normalize(ta - ro);
	float3 uu = normalize(cross(ww, float3(sin(roll), cos(roll), 0.0)));
	float3 vv = normalize(cross(uu, ww));
	return float3x3(uu, vv, ww);
}

float lerp(float a, float b, float x)
{
	return a + x * (b - a);
}

float3 lerp(float3 a, float3 b, float x)
{
	return float3(lerp(a.x, b.x, x), lerp(a.y, b.y, x), lerp(a.z, b.z, x));
}

float2 FoldX(float2 p)
{
	p.x = abs(p.x);
	return p;
}

float3 FoldX(float3 p)
{
	p.x = abs(p.x);
	return p;
}

float2x2 Rotate(float a)
{
	float s = sin(a), c = cos(a);
	return float2x2(c, s, -s, c);
}

float3x3 Rotate3D(float3 axis, float angle)
{
	//return float3x3(cos(angle) + axis.x * axis.x * (1.0 - cos(angle)), axis.x * axis.y * (1.0 - cos(angle))

	//手抜き
	return float3x3(1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0) * cos(angle)
		+ float3x3(axis.x * axis.x, axis.x * axis.y, axis.x * axis.z,
			axis.x * axis.y, axis.y * axis.y, axis.y * axis.z,
			axis.x * axis.z, axis.y * axis.z, axis.z * axis.z) * (1.0 - cos(angle))
		+ float3x3(0.0, -axis.z, axis.y,
			axis.z, 0.0, -axis.x,
			-axis.y, axis.x, 0.0) * sin(angle);
}

float3 TwistY(float3 p, float power)
{
	float s = sin(power * p.y);
	float c = cos(power * p.y);
	float3x3 m = float3x3(
		c, 0.0, -s,
		0.0, 1.0, 0.0,
		s, 0.0, c
	);

	return mul(m, p);
}

float2 iq_rand(float2 p)
{
	p = float2(dot(p, float2(127.1, 311.7)), dot(p, float2(269.5, 183.3)));
	return frac(sin(p)*43758.5453);
}

float3 nrand3(float2 co)
{
	float3 a = frac(cos(co.x*8.3e-3 + co.y)*float3(1.3e5, 4.7e5, 2.9e5));
	float3 b = frac(sin(co.x*0.3e-3 + co.y)*float3(8.1e5, 1.0e5, 0.1e5));
	float3 c = lerp(a, b, 0.5);
	return c;
}

float Sigmoid(float x, float gain, float offsetX)
{
	return (tanh(((x + offsetX)*gain) / 2.0) + 1.0) / 2.0;
}

//HSV -> RGB変換
//ネットに落ちてた実装そのまま
//https://qiita.com/masato_ka/items/c178a53c51364703d70b#_reference-6b1887cd1cd40aa64ace
float3 HSVtoRGB(float x)
{

	float gain = 10.0;
	float offsetX = 0.2;
	float offsetGreen = 0.6;

	x = (x * 2.0) - 1.0;
	float red = Sigmoid(x, gain, -1.0*offsetX);
	float blue = 1.0 - Sigmoid(x, gain, offsetX);
	float green = Sigmoid(x, gain, offsetGreen) + (1.0 - Sigmoid(x, gain, -1.0*offsetGreen));
	green = green - 1.0;
	return float3(blue, green, red);
}

float SoftMin(float a, float b, float r)
{
	float e = max(r - abs(a - b), 0.0);
	return min(a, b) - e * e*0.25 / r;
}

float SoftMax(float a, float b, float r)
{
	float e = max(r - abs(a - b), 0.0);
	return max(a, b) + e * e*0.25 / r;
}

float SmoothMin(in float a, in float b, in float k)
{
	return -log(exp(-k * a) + exp(-k * b)) / k;
}

float2 FoldRotate(in float2 p, in float s)
{
	float a = PI / s - atan2(p.x, p.y);
	float n = PI2 / s;
	a = floor(a / n) * n;
	p = mul(Rotate(a), p);
	return p;
}

float3 Mod(float3 a, float3 b)
{
	return frac(abs(a / b)) * abs(b);
}

float3 Repeat(float3 pos, float3 span)
{
	return Mod(pos, span) - span * 0.5;
}

float Union(float sdf1, float sdf2)
{
	return min(sdf1, sdf2);
}

float Intersection(float sdf1, float sdf2)
{
	return max(sdf1, sdf2);
}

float Subtraction(float sdf1, float sdf2)
{
	return max(sdf1, -sdf2);
}

//https://www.shadertoy.com/view/3dXXDN
//#TODO よくわかってない
float lengthN(float2 p, float n)
{
	p = pow(abs(p), float2(n, n));
	return pow(p.x + p.y, 1.0 / n);
}

//https://www.shadertoy.com/view/3dXXDN
//#TODO よくわかってない
float3 torusKnot(float t)
{
	t *= 6.283;
	float3 p = 0.3*float3(cos(t*7.0), sin(t*7.0), 0);
	p.x += 1.5;
	p.xz = mul(p.xz, Rotate(t*2.0));
	return p;
}

//----------------------------------------------------------------------------------------------
//Signed Distance Field
//----------------------------------------------------------------------------------------------
float SdfSphere(float3 currentRayPosition, float radius, float3 pos)
{
	// Here we get the distance to the surface of the balloon
	float distanceToBalloon = length(currentRayPosition - pos);

	// finally we get the distance to the balloon surface
	// by substacting the balloon radius. This means that if
	// the distance to the balloon is less than the balloon radius
	// the value we get will be negative! giving us the 'Signed' in
	// Signed Distance Field!
	return distanceToBalloon - radius;
}

float SdfBox(float3 currentRayPosition, float3 pos, float3 size)
{
	// Here we get the 'adjusted ray position' which is just
	// writing the point of the ray as if the origin of the 
	// space was where the box was positioned, instead of
	// at 0,0,0 . AKA the difference between the vectors in
	// vector format.
	float3 adjustedRayPosition = currentRayPosition - pos;

	// finally we get the distance to the box surface.
	float3 distanceVec = abs(adjustedRayPosition) - size;
	float maxDistance = max(distanceVec.x, max(distanceVec.y, distanceVec.z));
	float dist = min(maxDistance, 0.0); //箱の中にカメラがある場合？たぶん
	float distanceToBoxSurface = dist + length(max(distanceVec, 0.0));
	return distanceToBoxSurface;
}

float SdfTorus(float3 pos, float2 radius)
{
	float2 r = float2(length(pos.xy) - radius.x, pos.z);
	return length(r) - radius.y;
}

float SdfTorusKnot(float3 p)
{
	float ITR = 40.0, pitch = 1.0, t = 0.5, de = 1e10;
	for (int j = 0; j<2; j++)
	{
		float t0 = t - pitch * 0.5;
		pitch /= ITR;
		for (float i = 0.0; i <= ITR; i++)
		{
			t0 += pitch;
			float de0 = distance(p, torusKnot(t0));
			if (de0<de)
			{
				de = de0;
				t = t0;
			}
		}
	}

	float3 u = normalize(torusKnot(t));
	float3 v = normalize(torusKnot(t + 0.01) - torusKnot(t - 0.01));
	float3 w = normalize(cross(u, v));
	u = cross(v, w);
	p -= torusKnot(t);
	p = float3(dot(p, w), dot(p, u), dot(p, v));
	return lengthN(float2(length(p.yz), p.x), 3.0) - 0.18;
}

float SdfMandelbulb(in float3 p, in bool conservative)
{
    float3 w = p;
    float m = dot(w, w);

    float dz = 1.0;
    float3 J = float3(.2, 0.2, 0.2);
    
    float DERIVATIVE_BIAS = 1.0;

    for (int i = 0; i < 5; i++)
    {
        if (conservative)
            dz = max(dz * DERIVATIVE_BIAS, 8.0 * pow(m, 3.5) * dz + 1.0);
        else
            dz = 8.0 * pow(m, 3.5) * dz + 1.0;
		
        float r = length(w);
        float b = 8.0 * acos(clamp(w.y / r, -1.0, 1.0));
        float a = 8.0 * atan2(w.x, w.z);
        w = p + J + pow(r, 8.0) * float3(sin(b) * sin(a), cos(b), sin(b) * cos(a));
        
        m = dot(w, w);
		
        if (m > 4.0)
            break;
    }
       
    return 0.25 * log(m) * sqrt(m) / dz;
}

//-------------------------------------------------------------------------------------
// 'TAG : WHICH AM I CLOSER TO?'
// This function takes in two things
// and says which is closer by using the 
// distance to each thing, comparing them
// and returning the one that is closer!
float2 WhichThingAmICloserTo(float2 thing1, float2 thing2) {

	float2 closestThing = thing1;

	// Check out the balloon function
	// and remember how the x of the returned
	// information is the distance, and the y 
	// is the id of the thing!
	if (thing1.x <= thing2.x) {

		closestThing = thing1;

	}
	else if (thing2.x < thing1.x) {

		closestThing = thing2;

	}

	return closestThing;

}

#define MATERIAL_OPACITY 1.0
#define MATERIAL_PERFECT_REFLECTION 3.0
#define MATERIAL_ICE 10.0
#define MATERIAL_TEST_DIFFUSE 20.0
#define MATERIAL_TEST_DIFFUSE_LAMBERT 21.0
#define MATERIAL_TEST_SPECULAR_GGX 22.0
#define MATERIAL_TEST_REFLACT 23.0
#define MATERIAL_TEST_CAUSTICS 24.0
#define MATERIAL_TEST_CAUSTICS_GPUGEMS 25.0
#define MATERIAL_TEST_FAR 26.0

// Takes in the position of the ray, and feeds back
// 2 values of how close it is to things in the world
// what thing it is closest two in the world.
HitInfo MapTheWorld(in Ray ray, in float3 currentRayPosition)
{
	//Sphere
	float radius = 0.9;
	float3 spherePos = float3(0.0, 0.0, 0.0);
	float sdfSphere = SdfSphere(currentRayPosition, radius, spherePos);
	float2 sphere = float2(sdfSphere, MATERIAL_OPACITY);

	//Box
	float3 boxPosition = float3(0.0, 1.0, 0.0);
	float3 boxSize = float3(0.5, 0.5, 0.5);
	float sdfBox = SdfBox(mul(Rotate3D(YAxis, radians(45.0)), currentRayPosition), boxPosition, boxSize);
	//float sdfBox = SdfRoundBox(Rotate3D(YAxis, radians(45.0)) * currentRayPosition, boxPosition, boxSize, 0.01);
	float2 box = float2(sdfBox, MATERIAL_OPACITY);

	//TorusKnot
	float3 torusPos = float3(0.0, 0.0, 0.0);
	float2 torusRadius = float2(1.5, 0.12);
	float sdfTorus = SdfTorus(mul(Rotate3D(XAxis, radians(90.0)), currentRayPosition), torusRadius);
	float sdfTorusNot = SdfTorusKnot(currentRayPosition);
	float2 torusKnot = float2(min(sdfTorusNot, sdfTorus), MATERIAL_OPACITY);
	
	//MandelBulb
    float sdfMandelbulb = SdfMandelbulb(currentRayPosition, false);
    float2 mandelbulb = float2(sdfMandelbulb, MATERIAL_OPACITY);

	////Height Map Test
	//float2 uv = mod(abs(currentRayPosition.xz) * 0.25, 1.0);
	//float4 hightMap = texture(iChannel0, uv);
	//float heightFactor = 1.0;
	//hightMap *= heightFactor;

	//float3 normal = Rotate3D(ZAxis, radians(0.0)) * YAxis;
	//float sdfHeightMapPlane = SdfPlane(currentRayPosition, normal, 2.0);
	//float2 bumpedPlane = float2((sdfHeightMapPlane - hightMap.x * 1.0) * 0.1250, MATERIAL_OPACITY);

	////Water
	////float d = WaterSurface(currentRayPosition);   
	//float waveDistance = SdfPlane(currentRayPosition, YAxis, 2.0) - WaterWave(currentRayPosition);
	//float2 water = float2(waveDistance, MATERIAL_TEST_REFLACT);

	////Water With Caustics
	//float3 waterPos = currentRayPosition - float3(0.0, 0.75, 0.0);
	//float3 size = float3(2, 0, 2);
	//float waterDistance = length(max(abs(waterPos + float3(0, SdfPlane(waterPos, YAxis, 2.0) - WaterWave(waterPos), 0)) - size, 0.));
	//float2 waterCaustics = float2(waterDistance, MATERIAL_TEST_REFLACT);

	////Bottom With Caustics
	//float sdfBottomPlane = SdfPlane(currentRayPosition, YAxis, 2.0);
	//float2 bottom = float2(sdfBottomPlane, MATERIAL_TEST_CAUSTICS);

	////Bottom GPUGems
	//float2 bottom2 = float2(sdfBottomPlane, MATERIAL_TEST_CAUSTICS_GPUGEMS);

	////far
	//float farBallRadius = 1.0;
	//float3 farBallPos = float3(0.0, 0.0, 0.0);
	//float sdfFarBall = SdfSphere(currentRayPosition, farBallRadius, farBallPos);
	//float2 far = float2(sdfFarBall, MATERIAL_TEST_FAR);



    float2 result = WhichThingAmICloserTo(mandelbulb, mandelbulb);





	HitInfo info;
	info.distance = result.x;
	info.material = result.y;

	return info;




	//Fold Test
	//float3 treePos = currentRayPosition + float3(0.0, -0.0, -1.0);
	//float2 tree = float2(SdfSnowCrystal(treePos), 2);

	//Repeat Test
	//float2 repoatBox = SdfRepeatRoundBox(currentRayPosition);

	//Smin Test
	//float scale = 0.8;
	//float3 size = float3(0.1, 0.5, 0.1);
	//float3 smoothTreePos = currentRayPosition + float3(0.0, 1.5, 1.0);
	//float2 smoothTree = float2(SdfSmoothTree(smoothTreePos, scale, size, 8), 2);

	//Smooth Add Test
	//float3 posTorus = float3(0, 0.9, 0);
	//float2 radiusTorus = float2(0.75, 0.25);
	//float sdfSmoothTorusAndFloor = SmoothMin(SdfFloor(currentRayPosition), 
	//                                          SdfTorus((currentRayPosition - posTorus), radiusTorus),
	//                                          1.0);
	//float2 smoothTorusAndFloor = float2(sdfSmoothTorusAndFloor, 1.0);

	//Twist Test
	//float2 twistedYTorusRadius = float2(2.0, 0.6);
	//float3 twistedYTourusPos = currentRayPosition + float3(0.0, -0.0, 4.0);
	//float twist = 2.0; //g_time
	//float sdfTwistedYTorus = SdfTorus(TwistY(twistedYTourusPos, twist), twistedYTorusRadius);
	//float2 twistedYTorus = float2(sdfTwistedYTorus, 1.0);

	//RecursiveTetrahedron Test
	//float sdfRecursiveTetrahedron = SdfRecursiveTetrahedron(currentRayPosition);
	//float2 recursiveTetrahedron = float2(sdfRecursiveTetrahedron, 1.0);

	//Constructive Solid Geometry test
	//float3 csgBoxPos = float3( 0.0, -1.0, 0.0 );
	//float3 csgSize = float3(0.4, 0.4, 0.4);

	//float csgSphereRadius = .5;
	//float3 csgSpherePos = float3(0.0, -1.0, 0.0);
	//float sdfBoxCSG = SdfBox(currentRayPosition, csgBoxPos, csgSize);
	//float sdfShepereCSG = SdfSphere(currentRayPosition, csgSphereRadius, csgSpherePos);
	//float sdfSphereAndBoxUnion = Union(sdfBoxCSG, sdfShepereCSG);
	//float sdfSphereAndBoxIntersection = Intersection(sdfBoxCSG, sdfShepereCSG);    
	//float sdfSphereAndBoxSubtraction1 = Subtraction(sdfBoxCSG, sdfShepereCSG);
	//float sdfSphereAndBoxSubtraction2 = Subtraction(sdfShepereCSG, sdfBoxCSG);

	//float2 csgTest = float2(sdfSphereAndBoxIntersection, 2.0);

	//BarScene Test
	//float3 rotatedRay = currentRayPosition * Rotate3D(XAxis, radians(-45.0)) * Rotate3D(YAxis, radians(45.0));

	//float barWidth = 0.02;
	//float3 barInterval = float3(0.5);
	//float3 barRepeatPos = Repeat(rotatedRay, barInterval);
	//float sdfBarX = SdfBar(barRepeatPos.yz, barWidth);
	//float sdfBarY = SdfBar(barRepeatPos.xz, barWidth);
	//float sdfBarZ = SdfBar(barRepeatPos.xy, barWidth);    
	//float sdfBar = min(min(sdfBarX, sdfBarY), sdfBarZ);

	//float tubeWidth = 0.005;
	//float3 tubeInterval = float3(0.025);
	//float3 tubeRepeatPos = Repeat(rotatedRay, tubeInterval);
	//float sdfTubeX = SdfTube(tubeRepeatPos.yz, tubeWidth);
	//float sdfTubeY = SdfTube(tubeRepeatPos.xz, tubeWidth);
	//float sdfTubeZ = SdfTube(tubeRepeatPos.xy, tubeWidth);
	//float sdfTube = min(min(sdfTubeX, sdfTubeY), sdfTubeZ);

	//float2 barScene = float2(Subtraction(sdfBar, sdfTube), 2.0);

	//HexAnim Test　めんどくさくなってやめた
	//float2 hexSize = float2(1.0, 1.0);
	//float sdfHex = SdfHex(currentRayPosition.xz, hexSize);
	//float2 hex = float2(sdfHex, 2.0);

	//Morphing Test
	//float sdfMorphTest = SdfMorphSphereToRoundBox(currentRayPosition);
	//float2 morphTest = float2(sdfMorphTest, MATERIAL_OPACITY);

	//Advanced Animation
	//float sdfAnimeTest = SdfAdvancedBoxAnimation(currentRayPosition);
	//float2 animeTest = float2(sdfAnimeTest, MATERIAL_OPACITY);

	//Metaball(Sphere Tracing), SoftMin Test
	//float metaBallRadius = 0.5;
	//float3 metaBallPos1 = float3(.0, sin(g_time), 0.0);
	//float3 metaBallPos2 = float3(.0, cos(g_time), 0.0);
	//float metaBallSoftness = 0.5;
	//float metaBallPositiveSdf = SoftMin(SdfSphere(currentRayPosition, metaBallRadius, metaBallPos1),
	//                                      SdfSphere(currentRayPosition, metaBallRadius, metaBallPos2), metaBallSoftness);

	//float metaBallNegativeSdf = SoftMax(SdfSphere(currentRayPosition, metaBallRadius, metaBallPos1),
	//                                      -SdfSphere(currentRayPosition, metaBallRadius, metaBallPos2), metaBallSoftness);

	//float2 metaBall = float2(metaBallPositiveSdf, MATERIAL_TEST_REFLACT);

	//Terrain
	//float2 terrainUV = abs(currentRayPosition.xz + float2(9.0 + g_time, 10.0)); //10.0にすると真っ黒になる　よくわからん
	//float sdfTerrain = currentRayPosition.y + 1.0 - SdfTerrain(terrainUV) * 0.2;
	//float2 terrain = float2(sdfTerrain, MATERIAL_OPACITY);

}







//-------------------------------------------------------------------------------------
static const float HOW_CLOSE_IS_CLOSE_ENOUGH = 0.001;

// This is basically how big our scene is. each ray will be shot forward
// until it reaches this distance. the smaller it is, the quicker the 
// ray will reach the edge, which should help speed up this function
static const float FURTHEST_OUR_RAY_CAN_REACH = 20.;

// This is how may steps our ray can take. Hopefully for this
// simple of a world, it will very quickly get to the 'close enough' value
// and stop the iteration, but for more complex scenes, this value
// will dramatically change not only how good the scene looks
// but how fast teh scene can render. 
static const int HOW_MANY_STEPS_CAN_OUR_RAY_TAKE = 100;

// Here we are calcuting the normal of the surface
// Although it looks like alot of code, it actually
// is just trying to do something very simple, which
// is to figure out in what direction the SDF is increasing.
// What is amazing, is that this value is the same thing 
// as telling you what direction the surface faces, AKA the
// normal of the surface. 
float3 GetNormalOfSurface(in Ray ray, in float3 hitPos, float epsilon)
{
	// TODO 微小量を動的に決める実験　いまいち
	//epsilon *= length(ray.start - positionOfHit);

	return normalize(float3(
		MapTheWorld(ray, hitPos + XAxis * epsilon).distance - MapTheWorld(ray, hitPos - XAxis * epsilon).distance,
		MapTheWorld(ray, hitPos + YAxis * epsilon).distance - MapTheWorld(ray, hitPos - YAxis * epsilon).distance,
		MapTheWorld(ray, hitPos + ZAxis * epsilon).distance - MapTheWorld(ray, hitPos - ZAxis * epsilon).distance
	));
}

//dFdx,dFdyで法線を計算
//精度は低くなるが、わりとしっかり動く
float3 GetNormalOfSurface_byDf(in float3 positionOfHit)
{
	float3 dx = ddx(positionOfHit);
	float3 dy = ddy(positionOfHit);
	return normalize(cross(dx, dy));
}

HitInfo CheckRayHit(in Ray ray)
{
	float distanceToSurface = HOW_CLOSE_IS_CLOSE_ENOUGH * 2.;
	float totalDistanceTraveledByRay = 0.;
	float finalDistanceTraveledByRay = -1.;
	float finalID = -1.;
	float3 rayPos = float3(0.0, 0.0, 0.0);

	for (int i = 0; i < HOW_MANY_STEPS_CAN_OUR_RAY_TAKE; i++) {

		// First off, stop the iteration, if we are close enough to the surface!
		// Where more accurate close, far away is rough
		if (distanceToSurface < HOW_CLOSE_IS_CLOSE_ENOUGH) {
			//if( distanceToSurface < HOW_CLOSE_IS_CLOSE_ENOUGH * totalDistanceTraveledByRay * 0.01 ){
			break;
		}

		// Second off, stop the iteration, if we have reached the end of our scene! 
		if (totalDistanceTraveledByRay > FURTHEST_OUR_RAY_CAN_REACH) {
			break;
		}

		float3 currentPositionOfRay = ray.start + ray.dir * totalDistanceTraveledByRay;

		// Distance to and ID of things in the world
		HitInfo distanceAndIDOfThingsInTheWorld = MapTheWorld(ray, currentPositionOfRay);

		// we get out the results from our mapping of the world
		// I am reassigning them for clarity
		float distanceToThingsInTheWorld = distanceAndIDOfThingsInTheWorld.distance;
		float idOfClosestThingInTheWorld = distanceAndIDOfThingsInTheWorld.material;

		distanceToSurface = distanceToThingsInTheWorld;
		finalID = idOfClosestThingInTheWorld;
		rayPos = currentPositionOfRay;

		// ATTENTION: THIS THING IS AWESOME!
		// This last little calculation is probably the coolest hack
		// of this entire tutorial. If we wanted too, we could basically 
		// step through the field at a constant amount, and at every step    
		// say 'am i there yet', than move forward a little bit, and
		// say 'am i there yet', than move forward a little bit, and
		// that would take FOREVER, and get really annoying.

		// Instead what we say is 'How far until we are there?'
		// and move forward by that amount. This means that if
		// we are really far away from everything, we can make large
		// movements towards the surface, and if we are closer
		// we can make more precise movements. making our marching functino
		// faster, and ideally more precise!!

		// WOW!

		totalDistanceTraveledByRay += distanceToThingsInTheWorld;
	}

	// if we hit something set the finalDirastnce traveled by
	// ray to that distance!
	if (totalDistanceTraveledByRay < FURTHEST_OUR_RAY_CAN_REACH) {
		finalDistanceTraveledByRay = totalDistanceTraveledByRay;
	}
	else
	{
		finalID = 0.0;
	}

	float epsilon = 0.001;

	HitInfo info;
	info.ray = ray;
	info.distance = finalDistanceTraveledByRay;
	info.pos = rayPos;
	info.normal = GetNormalOfSurface(ray, rayPos, epsilon);
	//info.normal = GetNormalOfSurface_byDf(rayPos);
	info.material = finalID;

	return info;
}

//----------------------------------------------------------------------------------------------------------
//屈折率 Index of Refraction
#define IOR_AIR   1.0
#define IOR_WATER 1.3333333333
#define IOR_GLASS 1.5

// doing our background color is easy enough,
// just make it pure black. like my soul.
float3 doBackgroundColor() {
	return float3(0.0, 0.0, 0.0);
}

float FresnelSchlick(float dotVH, float f0)
{
	return f0 + (1.0 - f0) * (1.0 - pow(dotVH, 5.0));
}

float3 TestDiffuse(float3 positionOfHit, float3 normalOfSurface)
{

	float3 sunPosition = float3(1., 4., 3.);

	// the direction of the light goes from the sun
	// to the position of the hit
	float3 lightDirection = sunPosition - positionOfHit;


	// Here we are 'normalizing' the light direction
	// because we don't care how long it is, we
	// only care what direction it is!
	lightDirection = normalize(lightDirection);


	// getting the value of how much the surface
	// faces the light direction
	float faceValue = dot(lightDirection, normalOfSurface);

	// if the face value is negative, just make it 0.
	// so it doesn't give back negative light values
	// cuz that doesn't really make sense...
	faceValue = max(0., faceValue);

	float3 balloonColor = float3(1., 0., 0.);

	// our final color is the balloon color multiplied
	// by how much the surface faces the light
	float3 color = balloonColor * faceValue;

	// add in a bit of ambient color
	// just so we don't get any pure black
	color += float3(.3, .1, .2);


	return color;
}

//https://blog.selfshadow.com/publications/s2013-shading-course/
float BRDFDiffuseLambert(float dotVH)
{
	float value = 1.0 - pow(1.0 - dotVH, 5.0);
	return value / PI;
}

//http://jultika.oulu.fi/files/nbnfioulu-201805101776.pdf
float BRDFDiffuseLambertUE4(float dotNL)
{
	return dotNL / PI;
}

float Pow5(float x)
{
	return (x * x) * (x * x) * x;
}

// GGX Distribution
// Ref: https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
float DTerm(float3 l, float HDotN, float alpha2)
{
	float HDotN2 = HDotN * HDotN;
	float x = (1.0 - (1.0 - alpha2) * HDotN2);
	float D = alpha2 / (PI * x * x);
	return D;
}

// Smith Joint Masking-Shadowing Function
// Ref: https://hal.inria.fr/hal-00942452v1/document
float GTerm(float LDotN, float VDotN, float alpha2)
{
	float tanThetaLN2 = 1.0 / (LDotN * LDotN) - 1.0;
	float tanThetaVN2 = 1.0 / (VDotN * VDotN) - 1.0;

	float lambdaL = 0.5 * sqrt(1.0 + alpha2 * tanThetaLN2) - 0.5;
	float lambdaV = 0.5 * sqrt(1.0 + alpha2 * tanThetaVN2) - 0.5;

	return 1.0 / (1.0 + lambdaL + lambdaV);
}

// Schlick approx
// Ref: https://en.wikipedia.org/wiki/Schlick's_approximation
float FTerm(float LDotH, float f0)
{
	return f0 + (1.0 - f0) * Pow5(1.0 - LDotH);
}

float TorranceSparrowBRDF(float3 v, float3 l, float3 n, float roughness)
{
	float3 h = normalize(l + v);

	float VDotN = dot(v, n);
	float LDotN = dot(l, n);
	float HDotN = dot(h, n);
	float LDotH = dot(l, h);

	float alpha = roughness * roughness;
	float alpha2 = alpha * alpha;

	float D = DTerm(l, HDotN, alpha2);
	float G = GTerm(LDotN, VDotN, alpha2);
	float F = FTerm(LDotH, 0.95);
	return (D * G * F) / (4.0 * abs(VDotN) * abs(LDotN));
}

float3 SpecularGGX(HitInfo hit, in Scene scene)
{
	//mouse pos to roughness
	float2 iMouse = float2(0.5, 0.5);
	float2 iResolution = float2(1.0, 1.0);
	float2 mouseUV = iMouse.xy / iResolution.xy;
	if (iMouse.x <= 0.0) mouseUV = float2(0.2, 0.8);

	//roughness
	//Texture Sample
	float u = atan2(hit.pos.z, hit.pos.x) / (2.0 * PI);
	float v = atan2(hit.pos.y, length(hit.pos.xz)) / PI - 0.5;
	float2 tile = float2(10.0, 10.0);
	//float roughness = length(texture(iChannel2, tile * float2(u, v)).rgb) * mouseUV.y;
	float roughness = mouseUV.y;

	//color
	//float4 color = texture(iChannel2, tile * float2(u, v));
	float4 color = float4(1.0, 1.0, 1.0, 0.0);

	//radiance
	float radiance = TorranceSparrowBRDF(-hit.ray.dir, scene.dirLight.direction, hit.normal, roughness) * dot(scene.dirLight.direction, hit.normal);
	return float3(radiance * color.rgb) * 10.0;
}

// This is where we decide
// what color the world will be!
float3 ColorTheWorld(HitInfo hitInfo, Ray ray, Scene scene)
{
	float3 color = doBackgroundColor();

	float3 halfVector = normalize(scene.dirLight.direction + -hitInfo.ray.dir);
	float dotVH = max(0.0, dot(scene.dirLight.direction, halfVector));

	//放射照度 <- 放射輝度 * コサイン項
	float dotNL = max(0.0, dot(hitInfo.normal, scene.dirLight.direction));
	float3 irradiance = scene.dirLight.color * dotNL;

	//punctual lightの強さを調整
	//irradiance *= PI * scene.dirLight.scale; //scaleはてきとー　これでいいのか
	irradiance *= scene.dirLight.scale; //specularでしか使わなくなったので消した

	//roughness
	//mouse pos to roughness
	float2 iMouse = float2(0.5, 0.5);
	float2 iResolution = float2(1.0, 1.0);
	float2 mouseUV = iMouse.xy / iResolution.xy;
	if (iMouse.x <= 0.0) mouseUV = float2(0.2, 0.8);

	if (hitInfo.material == MATERIAL_OPACITY) {
		//roughness from texture
		float u = atan2(hitInfo.pos.z, hitInfo.pos.x) / (2.0 * PI);
		float v = atan2(hitInfo.pos.y, length(hitInfo.pos.xz)) / PI - 0.5;
		float2 tile = float2(10.0, 10.0);
		//float roughness = length(texture(iChannel2, tile * float2(u, v)).rgb) * mouseUV.y;
		float roughness = mouseUV.y;

		//metallic
		float metallic = mouseUV.x;

		//color
		//float3 baseColor = texture(iChannel2, tile * float2(u, v)).rgb;
		float3 baseColor = float3(1, 1, 1);

		//Diffuse
		float diffuseBRDF = BRDFDiffuseLambertUE4(1.0) * dotNL;
		//float diffuseBRDF = BRDFDiffuseLambert(1.0) * dotNL;
		float3 diffuseColor = diffuseBRDF * baseColor;
		//float3 diffuseColor = baseColor * (1.0 - metallic) ;

		//Specular
		float specularBRDF = TorranceSparrowBRDF(-hitInfo.ray.dir, scene.dirLight.direction, hitInfo.normal, roughness);
		float3 specularColor = specularBRDF * irradiance * metallic;

		color = diffuseColor + specularColor;

	}
	else if (hitInfo.material == MATERIAL_PERFECT_REFLECTION) {
		//Perfect Reflection
		float3 reflected = reflect(ray.dir, hitInfo.normal);
		//color = texture(iChannel1, reflected).rgb;

	}
	else if (hitInfo.material == MATERIAL_ICE) {
		float refractionIndex = mouseUV.x + 1.0;

		//Reflect
		float3 reflected = reflect(ray.dir, hitInfo.normal);
		float3 refracted = refract(ray.dir, hitInfo.normal, 1.0 / refractionIndex);
		float lightDotReflect = max(0.0, dot(reflected, scene.dirLight.color));

		//二回屈折　ちょっとてきとー
		Ray reflectedRay = { hitInfo.pos, reflected };
		HitInfo secondHit = CheckRayHit(reflectedRay);
		reflected = reflect(secondHit.ray.dir, secondHit.normal);
		refracted = refract(secondHit.ray.dir, secondHit.normal, 1.0 / refractionIndex);

		float f0 = pow(1.0 - refractionIndex, 2.0) / (1.0 + refractionIndex);
		float fresnel = FresnelSchlick(dotVH, f0);
		float3 specular = fresnel * scene.dirLight.color * lightDotReflect;

		//てきとーな調整
		color = /*texture(iChannel1, refracted).rgb + */ float3(0.0, 0.0, 0.0) + specular * 0.2 + float3(0.0, 0.00, 0.06);

	}
	else if (hitInfo.material == MATERIAL_TEST_REFLACT) {
		float refractionIndex = mouseUV.x + 1.0;

		//Reflect
		float3 reflected = reflect(ray.dir, hitInfo.normal);
		float3 refracted = refract(ray.dir, hitInfo.normal, 1.0 / refractionIndex);
		float lightDotReflect = max(0.0, dot(reflected, scene.dirLight.color));

		Ray refractedRay = { hitInfo.pos, refracted };
		HitInfo secondHit = CheckRayHit(refractedRay);
		reflected = reflect(secondHit.ray.dir, secondHit.normal);
		refracted = refract(secondHit.ray.dir, secondHit.normal, 1.0 / refractionIndex);

		float f0 = pow(1.0 - refractionIndex, 2.0) / (1.0 + refractionIndex);
		float fresnel = FresnelSchlick(dotVH, f0);
		float3 specular = fresnel * scene.dirLight.color * lightDotReflect;

		//二回屈折　ちょっとてきとー
		//color = ColorTheWorld(secondHit, refractedRay, scene);
		color = /*texture(iChannel1, refracted).rgb*/float3(0.0, 0.0, 0.0) + specular;

	}
	else if (hitInfo.material == MATERIAL_TEST_DIFFUSE) {
		//Simple Diffuse        
		color = TestDiffuse(hitInfo.pos, hitInfo.normal);

	}
	else if (hitInfo.material == MATERIAL_TEST_DIFFUSE_LAMBERT) {
		//Lambert
		float3 balloonColor = float3(0.5, 0.5, 0.5);
		float3 testColor = balloonColor * BRDFDiffuseLambert(1.0);
		//punctual lightの強さを調整
		irradiance *= PI; //scaleはてきとー　これでいいのか
		color = irradiance * testColor;

	}
	else if (hitInfo.material == MATERIAL_TEST_SPECULAR_GGX) {
		//GGX
		color = SpecularGGX(hitInfo, scene);

	}
	else if (hitInfo.material == MATERIAL_TEST_CAUSTICS) {
		////Caustics
		////無理やり実装してみる
		////二つの平行な平面が上と下で並んでいる、とする
		////上は水面、下は底
		////底に当たる光を計算してみる


		////平面からカメラ方向に一度レイを飛ばして、水面の位置を計算

		////光源からライティングする水底へのベクトル　位置はてきとー
		//float3 lightPos = float3(0.0, 10.0, 0.0); //水面よりも遠くないと変なことになる
		//float3 lightRay = normalize(hitInfo.pos - lightPos);

		////変形前の水面の位置　位置はてきとー
		//float3 stopSurfacePos = hitInfo.pos - float3(0.0, 1.0, 0.0);

		////変形前の水面を屈折したレイ
		//float ior = IOR_AIR / IOR_WATER;
		//float3 stopSurfaceRefractedRay = refract(lightRay, YAxis, ior); //法線はてきとーに上方向

		//															  //屈折したレイと水底との交点
		//															  //TODO 実装間違ってる…　正しい底の位置が取れてない
		//float stopHit = SdfPlane(stopSurfacePos, YAxis, 2.0);
		//float3 beforePos = stopSurfacePos + stopSurfaceRefractedRay * stopHit;

		////変形後の水面の位置
		//float3 waveSurfacePos = stopSurfacePos + float3(0.0, WaterWave(stopSurfacePos), 0.0);

		////変形後の水面で光線を屈折
		//float epsilon = 0.0001;
		//Ray afterRay = Ray(lightPos, lightRay);
		//float3 waveNormal = GetNormalOfSurface(afterRay, waveSurfacePos, epsilon); //うまく計算できてないっぽい
		//float3 afterRefractedRay = refract(lightRay, waveNormal, ior);
		//float afterHit = SdfPlane(waveSurfacePos, YAxis, 2.0);
		//float3 afterPos = waveSurfacePos + afterRefractedRay * stopHit;

		////変形前後での面積の変化量を明るさの比率に
		//float beforeArea = length(dFdx(beforePos)) * length(dFdy(beforePos));
		//float afterArea = length(dFdx(afterPos)) * length(dFdy(afterPos));

		////color = float3(0.0, 0.0, 0.5) * max(beforeArea/afterArea, 0.001);
		//color = float3(max(beforeArea / afterArea, 0.001) * 0.25);

		//float2 uv = abs(hitInfo.pos.xz) * 0.125;
		//color += float3(texture(iChannel2, uv).rgb) * 0.5;
		////color = waveNormal;

	}
	else if (hitInfo.material == MATERIAL_TEST_CAUSTICS_GPUGEMS) {
		//// Generate a normal (line direction) from the gradient

		//// of the wave function and intercept with the water plane.

		//// We use screen-space z to attenuate the effect to avoid aliasing.

		//float3 wavePos = float3(0.0, 0.0, 0.0);
		//float distance = 0.0;
		//float2 dxdy = gradwave(wavePos.x, wavePos.y, g_time);

		//float3 intercept = line_plane_intercept(
		//	wavePos,
		//	float3(dxdy, clamp(0.0, 1.0, distance)),
		//	float3(0, 0, 1), -0.8);

		//// OUTPUT

		//float3 colour;
		//float2 uv = abs(hitInfo.pos.xz) * 0.5;
		//colour.rgb += float3(1.0 - length(intercept.xy)); //light map
		//												//colour.rgb += float3(texture(iChannel2, uv).rgb);
		//return colour;
	}
	else if (hitInfo.material == 0.0) {
		//Cubemap
		//color = texture(iChannel1, ray.dir).rgb;
		color = texCube.Sample(sampler0, ray.dir);
	}
	else if (hitInfo.material == MATERIAL_TEST_FAR) {
		//float3 p = float3(0.0, 0.0, 0.0);
		//const float r = 1.0;
		//float t;
		//float4 c = float4(0.0, 0.0, 0.0, 0.0);
		//float3 pos = hitInfo.pos + ray.dir * 0.075;

		//// ray-march into volume
		//for (int i = 0; i < furLayers; i++) {
		//	float4 sampleCol;
		//	float2 uv;
		//	sampleCol.a = furDensity(pos, uv);
		//	if (sampleCol.a > 0.0) {
		//		sampleCol.rgb = furShade(pos, uv, ray.start, sampleCol.a);

		//		// pre-multiply alpha
		//		sampleCol.rgb *= sampleCol.a;
		//		c = c + sampleCol * (1.0 - c.a);
		//		if (c.a > 0.95) break;
		//	}

		//	pos += ray.dir*rayStep;
		//}

		//color = c.rgb;
	}

	return color;
}


//-------------------------------------------------------------------------------------
float3x3 CalculateEyeRayTransformationMatrix(in float3 ro, in float3 ta, in float roll)
{
	float3 ww = normalize(ta - ro);
	float3 uu = normalize(cross(ww, float3(sin(roll), cos(roll), 0.0)));
	float3 vv = normalize(cross(uu, ww));
	return float3x3(uu, vv, ww);
}

float2 NormalizingCoordinated(in float2 fragCoord, in float2 resolution)
{
	// Here we are getting our 'Position' of each pixel
	// This section is important, because if we didn't
	// divied by the resolution, our values would be masssive
	// as fragCoord returns the value of how many pixels over we are. which is alot :)

	//[-1,1] along the shortest side
	return (2.0 * fragCoord.xy - resolution.xy) / min(resolution.x, resolution.y);

	//[0,1]
	//return fragCoord / min(iResolution.x, iResolution.y);    
}

//-------------------------------------------------------------------------------------
PSInput VSMain(VSInput input)
{
	PSInput	result;

	result.position = input.pos;
	result.uv = float2(input.uv.x, 1.0f - input.uv.y);	//ShaderToyの座標系に合わせた　もちろんcppの頂点側の初期値で調整してもよい

	return result;
}

float4 PSMain(PSInput input) : SV_TARGET
{
	//#define WINDOW_WIDTH	1280
	//#define WINDOW_HEIGHT	720
	//float	g_aspectRatio = (float)WINDOW_WIDTH / (float)WINDOW_HEIGHT;
	//float2 inScreenPos = NormalizingCoordinated(input.uv, float2(1, 1/g_aspectRatio)) - float2(0.0, 0.5);
	
	float2 inScreenPos = NormalizingCoordinated(input.uv, float2(1, 1));

	//move camera
	//float angle = 0.2 * g_time; float3 camPos = float3(10.0 * sin(angle), 2.5 * cos(0.4 * angle), 10.0 * cos(angle));

	//rotate camera
	//float3 camPos = float3(1.5 * sin(1.5 * g_time), 1.0, 6.0);

	float3 camPos = float3(0., 1.0, 2.);
	float3 camLookAt = float3(0., 0, 0.);
	float3x3 eyeTransformationMtx = CalculateEyeRayTransformationMatrix(camPos, camLookAt, 0.);

	//scene
	//DirectionalLight dirLight = DirectionalLight(normalize(float3(1., 1., 1.)), float3(1.0, 1.0, 1.0), 1.0);
	//Scene scene = Scene(dirLight);
	DirectionalLight dirLight = { normalize(float3(1., 1., 1.)), float3(1.0, 1.0, 1.0), 1.0 };
	Scene scene = { dirLight };

	// Here we get the actual ray that goes out of the eye
	// and through the individual pixel! This basically the only thing
	// that is different between the pixels, but is also the bread and butter of ray tracing. 
	// It should be since it has the word 'ray' in its variable name...
	// the 2. at the end is the 'lens length' . I don't know how to best describe this,
	// but once the full scene is built, tryin playing with it to understand inherently how it works
	float3 rayComingOutOfEyeDirection = normalize(mul(eyeTransformationMtx, float3(inScreenPos.xy, 2.0f)));

	//通常描画
	Ray ray = { camPos, rayComingOutOfEyeDirection };
	HitInfo hitInfo = CheckRayHit(ray);
	float3 color = ColorTheWorld(hitInfo, ray, scene);

	//MSAAもどき　4サンプル
	float epsilon = 0.001;
	//color = MSAA(scene, ray, color, hitInfo, epsilon);

	//color = float4(DebugRenderMode(color, hitInfo, color, DEBUG_RENDER_MODE_NONE), 1.0);
	return float4(color, 0.0);
	//return tex0.Sample(sampler0, input.uv);
	//return tex1.Sample(sampler0, input.uv);
	//return texCube.Sample(sampler0, ray.dir);

	//#TODO デバッグモードとしてうまく扱いたい
	//return float4(inScreenPos, 0.0, 0.0);
	//return float4(input.uv, 0.0, 0.0);
}