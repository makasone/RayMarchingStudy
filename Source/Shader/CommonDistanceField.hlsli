#include "CommonFunction.hlsli"


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
    float dist = min(maxDistance, 0.0); //î†ÇÃíÜÇ…ÉJÉÅÉâÇ™Ç†ÇÈèÍçáÅHÇΩÇ‘ÇÒ
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
    for (int j = 0; j < 2; j++)
    {
        float t0 = t - pitch * 0.5;
        pitch /= ITR;
        for (float i = 0.0; i <= ITR; i++)
        {
            t0 += pitch;
            float de0 = distance(p, torusKnot(t0));
            if (de0 < de)
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

float SdfHex(float2 pos, float2 height)
{
    float2 q = abs(pos);
    return max(q.x + q.y * 0.57735, q.y * 1.1547) - height.x;
}