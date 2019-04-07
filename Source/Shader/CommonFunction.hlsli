//----------------------------------------------------------------------------------------------
#define PI	3.14159265359
#define PI2	PI * 2.0
#define XAxis float3(1.0, 0.0, 0.0)
#define YAxis float3(0.0, 1.0, 0.0)
#define ZAxis float3(0.0, 0.0, 1.0)
//----------------------------------------------------------------------------------------------


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

	//Žè”²‚«
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
    return frac(sin(p) * 43758.5453);
}

float3 nrand3(float2 co)
{
    float3 a = frac(cos(co.x * 8.3e-3 + co.y) * float3(1.3e5, 4.7e5, 2.9e5));
    float3 b = frac(sin(co.x * 0.3e-3 + co.y) * float3(8.1e5, 1.0e5, 0.1e5));
    float3 c = lerp(a, b, 0.5);
    return c;
}

float Sigmoid(float x, float gain, float offsetX)
{
    return (tanh(((x + offsetX) * gain) / 2.0) + 1.0) / 2.0;
}

//HSV -> RGB•ÏŠ·
//ƒlƒbƒg‚É—Ž‚¿‚Ä‚½ŽÀ‘•‚»‚Ì‚Ü‚Ü
//https://qiita.com/masato_ka/items/c178a53c51364703d70b#_reference-6b1887cd1cd40aa64ace
float3 HSVtoRGB(float x)
{

    float gain = 10.0;
    float offsetX = 0.2;
    float offsetGreen = 0.6;

    x = (x * 2.0) - 1.0;
    float red = Sigmoid(x, gain, -1.0 * offsetX);
    float blue = 1.0 - Sigmoid(x, gain, offsetX);
    float green = Sigmoid(x, gain, offsetGreen) + (1.0 - Sigmoid(x, gain, -1.0 * offsetGreen));
    green = green - 1.0;
    return float3(blue, green, red);
}

float SoftMin(float a, float b, float r)
{
    float e = max(r - abs(a - b), 0.0);
    return min(a, b) - e * e * 0.25 / r;
}

float SoftMax(float a, float b, float r)
{
    float e = max(r - abs(a - b), 0.0);
    return max(a, b) + e * e * 0.25 / r;
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
//#TODO ‚æ‚­‚í‚©‚Á‚Ä‚È‚¢
float lengthN(float2 p, float n)
{
    p = pow(abs(p), float2(n, n));
    return pow(p.x + p.y, 1.0 / n);
}

//https://www.shadertoy.com/view/3dXXDN
//#TODO ‚æ‚­‚í‚©‚Á‚Ä‚È‚¢
float3 torusKnot(float t)
{
    t *= 6.283;
    float3 p = 0.3 * float3(cos(t * 7.0), sin(t * 7.0), 0);
    p.x += 1.5;
    p.xz = mul(p.xz, Rotate(t * 2.0));
    return p;
}


// A small tweak to iq's mandelbulb distance estimator
// that considers Julia offsets. Right side is better.

// The main problem arises when the floating bulbs appear, 
// as there's an invisible sphere of overstepping that prevents rendering the outer bulbs.

// The tweak is not perfect, and you can force more detail by increasing DERIVATIVE_BIAS to > 1.0
// Beware that iteration count will increase.

// In the left side, you should see different banding patterns, and some far bulbs disappearing.
// In the right side, the banding disappears and bulbs are better estimated.

// Note that this is not just biasing the result, because this "sphere of death" is also fractal-like,
// and by modifying the derivative it properly adapts to the solution.

// Color is iteration count; red is more.

// I'll add more later, as I currently lack the time and coffee to rationalize this result.

// More info on http://www.fractalforums.com/new-theories-and-research/error-estimation-of-distance-estimators/

#define MAX_ITERATIONS 200
#define EPSILON .0001
#define DERIVATIVE_BIAS 1.0

#define saturate(x) clamp(x, 0.0, 1.0)

float evaluateMandelbulb(in vec3 p, in bool conservative)
{
    vec3 w = p;
    float m = dot(w, w);

    float dz = 1.0;
    vec3 J = vec3(.2);
       
    for (int i = 0; i < 5; i++)
    {
        if (conservative)
            dz = max(dz * DERIVATIVE_BIAS, 8.0 * pow(m, 3.5) * dz + 1.0);
        else
            dz = 8.0 * pow(m, 3.5) * dz + 1.0;
		
        float r = length(w);
        float b = 8.0 * acos(clamp(w.y / r, -1.0, 1.0));
        float a = 8.0 * atan(w.x, w.z);
        w = p + J + pow(r, 8.0) * vec3(sin(b) * sin(a), cos(b), sin(b) * cos(a));
        
        m = dot(w, w);
		
        if (m > 4.0)
            break;
    }
       
    return 0.25 * log(m) * sqrt(m) / dz;
}

float trace(vec3 rayO, vec3 rayD, bool conservative)
{
    float t = 0.0;
    float steps = 0.0;
    float it = 1.0 / float(MAX_ITERATIONS);
    
    for (int j = 0; j < MAX_ITERATIONS; j++)
    {
        vec3 p = rayO + rayD * t;
        float d = evaluateMandelbulb(p, conservative);

        if (d < EPSILON)
            break;

        t += d;
        
        if (t >= .5)
            break;
        
        steps += it;
    }
    
    return steps;
}

// Reference: http://www.iquilezles.org/www/articles/palettes/palettes.htm
vec3 palette(float t, vec3 a, vec3 b, vec3 c, vec3 d)
{
    return saturate(a + b * cos(6.28318 * ( c * t + d)));
}

vec3 debugIterations(float factor)
{
    vec3 a = vec3(0.478, 0.500, 0.500);
    vec3 b = vec3(0.500);
    vec3 c = vec3(0.688, 0.748, 0.748);
    vec3 d = vec3(0.318, 0.588, 0.908);

    return palette(factor, a, b, c, d);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
    
    float time = iTime * .2;
    vec3 target = vec3(.9, 0.22, -0.10);
    vec3 p = target + vec3(cos(time), 0.0, sin(time)) * .15;
    
    vec3 forward = normalize(target - p);
    vec3 left = normalize(cross(forward, vec3(0.0, 1.0, 0.0)));
    vec3 up = normalize(cross(forward, left));
    
    vec3 rayOrigin = p;
    vec3 rayDirection = normalize(forward + left * uv.x - up * uv.y);
    
    bool conservative = false;
    float threshold = (-iResolution.x + 2.0 * iMouse.x) / iResolution.y;
    
    if (iMouse.z < 0.01)
        threshold = 0.5 / iResolution.y;
    
    if (uv.x > threshold)
        conservative = true;

    float steps = trace(rayOrigin, rayDirection, conservative);
    
    vec3 color = debugIterations(steps);
    color = mix(vec3(0.0), color, smoothstep(0.002, .005, abs(uv.x - threshold)));
    
    fragColor = vec4(color, 1.0);
}