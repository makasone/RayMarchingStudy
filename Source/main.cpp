#include <windows.h>
#include <tchar.h>
#include <wrl.h>		// Microsoft::WRL::ComPtr
#include <vector>

#include <d3d12.h>
#pragma comment(lib, "d3d12.lib")
#include <dxgi1_4.h>
#pragma comment(lib, "dxgi.lib")
#include <D3Dcompiler.h>
#pragma comment(lib, "d3dcompiler.lib")
#include <DirectXMath.h>

#include "TexReader/DirectXTex.h"

using namespace DirectX;
using namespace Microsoft::WRL;

#define WINDOW_CLASS	_T("DirectX12Test")
#define WINDOW_TITLE	WINDOW_CLASS
#define	WINDOW_STYLE	WS_OVERLAPPEDWINDOW
//#define WINDOW_WIDTH	1280
//#define WINDOW_HEIGHT	720
#define WINDOW_WIDTH	512
#define WINDOW_HEIGHT	512

// ウィンドウプロシージャ
LRESULT CALLBACK WindowProc(HWND hWnd, UINT nMsg, WPARAM wParam, LPARAM lParam);
// 初期化
BOOL Init(HWND hWnd);
BOOL InitTexture(HWND hWnd);
// 描画
BOOL Draw();
// 
BOOL WaitForPreviousFrame();

const UINT	FrameCount = 2;

struct Vertex {
	XMFLOAT3	position;
	XMFLOAT2	uv;
};

//--------------------------------------------------------------------------------------
#pragma pack(push,1)
struct SimpleVertex
{
	XMFLOAT4 Pos;
	XMFLOAT4 Tex;
};

struct CBArrayControl
{
	float Index;
	float pad[3];
};
#pragma pack(pop)

D3D_FEATURE_LEVEL					g_featureLevel = D3D_FEATURE_LEVEL_11_0;

// パイプラインオブジェクト
D3D12_VIEWPORT						g_viewport = { 0.0f, 0.0f, (float)WINDOW_WIDTH, (float)WINDOW_HEIGHT, 0.0f, 1.0f };
D3D12_RECT							g_scissorRect = { 0, 0, WINDOW_WIDTH, WINDOW_HEIGHT };
ComPtr<IDXGISwapChain3>				g_swapChain;
ComPtr<ID3D12Device>				g_device;
ComPtr<ID3D12Resource>				g_renderTargets[FrameCount];
ComPtr<ID3D12CommandAllocator>		g_commandAllocator;
ComPtr<ID3D12CommandQueue>			g_commandQueue;
ComPtr<ID3D12RootSignature>			g_rootSignature;
ComPtr<ID3D12DescriptorHeap>		g_descriptorHeap;
ComPtr<ID3D12PipelineState>			g_pipelineState;
ComPtr<ID3D12GraphicsCommandList>	g_commandList;
UINT								g_descriptorSize = 0;
UINT								g_cbvSrvDescriptorSize = 0;

//リソース
ComPtr<ID3D12Resource>		g_vertexBuffer;
D3D12_VERTEX_BUFFER_VIEW	g_vertexBufferView;
ComPtr<ID3D12Resource>		g_indexBuffer;
D3D12_INDEX_BUFFER_VIEW		g_indexBufferView;

//テクスチャー
ComPtr<ID3D12Resource>			g_rTexture;
ComPtr<ID3D12Resource>			g_rTestTexture;
ComPtr<ID3D12Resource>			g_rCubeMapTexture;
ComPtr<ID3D12DescriptorHeap>	g_dhTexture;
ComPtr<ID3D12GraphicsCommandList>	g_texCommandList;

//定数バッファ
ComPtr<ID3D12Resource>			g_constantBuffer;

// 同期オブジェクト
UINT				g_frameIndex = 0;
HANDLE				g_fenceEvent;
ComPtr<ID3D12Fence>	g_fence;
UINT64				g_fenceValue;

// ビューポートのアスペクト比
float	g_aspectRatio = (float)WINDOW_WIDTH / (float)WINDOW_HEIGHT;

//シーン情報
float g_time = 0.0f;

// アダプタ情報
bool	g_useWarpDevice = false;
//bool	g_useWarpDevice = true;


int WINAPI _tWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, TCHAR *lpszCmdLine, int nCmdShow)
{
	// ウィンドウを作成
	WNDCLASSEX	wndclass = {};

	// ウィンドウクラスを登録
	wndclass.cbSize = sizeof(WNDCLASSEX);
	wndclass.style = CS_HREDRAW | CS_VREDRAW | CS_DBLCLKS;
	wndclass.lpfnWndProc = WindowProc;
	wndclass.hInstance = hInstance;
	wndclass.hCursor = LoadCursor(NULL, IDC_ARROW);
	wndclass.lpszClassName = WINDOW_CLASS;
	RegisterClassEx(&wndclass);

	RECT	windowRect = { 0, 0, WINDOW_WIDTH, WINDOW_HEIGHT };

	AdjustWindowRect(&windowRect, WINDOW_STYLE, FALSE);

	// ウィンドウを作成
	HWND	hWnd = CreateWindow(
		WINDOW_CLASS,
		WINDOW_TITLE,
		WINDOW_STYLE,
		CW_USEDEFAULT,
		CW_USEDEFAULT,
		windowRect.right - windowRect.left,
		windowRect.bottom - windowRect.top,
		NULL,
		NULL,
		hInstance,
		NULL);

	// DirectXを初期化
	if (!Init(hWnd)) {
		MessageBox(hWnd, _T("DirectXの初期化が失敗しました"), _T("Init"), MB_OK | MB_ICONEXCLAMATION);
		return 0;
	}

	ShowWindow(hWnd, SW_SHOW);
	UpdateWindow(hWnd);

	// メッセージループ
	MSG	msg;

	while (1) {
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
			if (msg.message == WM_QUIT) break;

			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}

	// 終了時の後処理
	WaitForPreviousFrame();
	CloseHandle(g_fenceEvent);

	return (int)msg.wParam;
}

// ウィンドウプロシージャ
LRESULT CALLBACK WindowProc(HWND hWnd, UINT nMsg, WPARAM wParam, LPARAM lParam)
{
	switch (nMsg) {
	case WM_PAINT:
		// 描画
		Draw();
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hWnd, nMsg, wParam, lParam);
	}

	return 0;
}

// 初期化
BOOL Init(HWND hWnd)
{
#if defined(_DEBUG)
	// DirectX12のデバッグレイヤーを有効にする
	{
		ComPtr<ID3D12Debug>	debugController;
		if (SUCCEEDED(D3D12GetDebugInterface(IID_PPV_ARGS(&debugController)))) {
			debugController->EnableDebugLayer();
		}
	}
#endif

	// DirectX12がサポートする利用可能なハードウェアアダプタを取得
	ComPtr<IDXGIFactory4>	factory;
	if (FAILED(CreateDXGIFactory1(IID_PPV_ARGS(&factory)))) return FALSE;

	//#TODO デフォルトがなかったらfall backしたほうがよい？
	if (g_useWarpDevice) {
		ComPtr<IDXGIAdapter>	warpAdapter;
		if (FAILED(factory->EnumWarpAdapter(IID_PPV_ARGS(&warpAdapter)))) return FALSE;
		if (FAILED(D3D12CreateDevice(warpAdapter.Get(), g_featureLevel, IID_PPV_ARGS(&g_device)))) return FALSE;
	}
	else {
		ComPtr<IDXGIAdapter1>	hardwareAdapter;
		ComPtr<IDXGIAdapter1>	adapter;
		hardwareAdapter = nullptr;

		for (UINT i = 0; DXGI_ERROR_NOT_FOUND != factory->EnumAdapters1(i, &adapter); i++) {
			DXGI_ADAPTER_DESC1 desc;
			adapter->GetDesc1(&desc);
			if (desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE) continue;
			// アダプタがDirectX12に対応しているか確認
			if (SUCCEEDED(D3D12CreateDevice(adapter.Get(), g_featureLevel, _uuidof(ID3D12Device), nullptr))) break;
		}

		hardwareAdapter = adapter.Detach();

		if (FAILED(D3D12CreateDevice(hardwareAdapter.Get(), g_featureLevel, IID_PPV_ARGS(&g_device)))) return FALSE;
	}

	// コマンドキューを作成
	D3D12_COMMAND_QUEUE_DESC	queueDesc = {};
	queueDesc.Flags = D3D12_COMMAND_QUEUE_FLAG_NONE;
	queueDesc.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;

	if (FAILED(g_device->CreateCommandQueue(&queueDesc, IID_PPV_ARGS(&g_commandQueue)))) return FALSE;

	// スワップチェインを作成
	DXGI_SWAP_CHAIN_DESC1	swapChainDesc = {};
	swapChainDesc.BufferCount = FrameCount;
	swapChainDesc.Width = WINDOW_WIDTH;
	swapChainDesc.Height = WINDOW_HEIGHT;
	swapChainDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
	swapChainDesc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	swapChainDesc.SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD;
	swapChainDesc.SampleDesc.Count = 1;

	ComPtr<IDXGISwapChain1>	swapChain;
	if (FAILED(factory->CreateSwapChainForHwnd(g_commandQueue.Get(), hWnd, &swapChainDesc, nullptr, nullptr, &swapChain))) return FALSE;

	// フルスクリーンのサポートなし
	if (FAILED(factory->MakeWindowAssociation(hWnd, DXGI_MWA_NO_ALT_ENTER))) return FALSE;

	if (FAILED(swapChain.As(&g_swapChain))) return FALSE;
	g_frameIndex = g_swapChain->GetCurrentBackBufferIndex();

	// 記述子ヒープを作成
	{
		// レンダーターゲットビュー用の記述子ヒープを作成
		D3D12_DESCRIPTOR_HEAP_DESC	descriptorHeapDesc = {};
		descriptorHeapDesc.NumDescriptors = FrameCount;
		descriptorHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;
		descriptorHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
		if (FAILED(g_device->CreateDescriptorHeap(&descriptorHeapDesc, IID_PPV_ARGS(&g_descriptorHeap)))) return FALSE;

		g_descriptorSize = g_device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);
	}

	// フレームリソースを作成
	{
		D3D12_CPU_DESCRIPTOR_HANDLE	cpuDescriptorHandle = {};
		cpuDescriptorHandle.ptr = g_descriptorHeap->GetCPUDescriptorHandleForHeapStart().ptr;

		// フレームバッファとバックバッファののレンダーターゲットビューを作成
		for (UINT i = 0; i < FrameCount; i++) {
			if(FAILED(g_swapChain->GetBuffer(i, IID_PPV_ARGS(&g_renderTargets[i])))) return FALSE;
			g_device->CreateRenderTargetView(g_renderTargets[i].Get(), nullptr, cpuDescriptorHandle);
			cpuDescriptorHandle.ptr += g_descriptorSize;
		}
	}

	// コマンドアロケーターを作成
	if (FAILED(g_device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&g_commandAllocator)))) return FALSE;

	// ルートシグネチャを作成
	{
		D3D12_DESCRIPTOR_RANGE	descriptorRange[1] = {};

		//テクスチャを使うための設定
		descriptorRange[0].NumDescriptors = 1 + 1 ;
		descriptorRange[0].BaseShaderRegister = 0;
		descriptorRange[0].RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_SRV;
		descriptorRange[0].OffsetInDescriptorsFromTableStart = D3D12_DESCRIPTOR_RANGE_OFFSET_APPEND;

		//パラメータ
		D3D12_ROOT_PARAMETER	rootParameters[2] = {};
		rootParameters[0].ParameterType = D3D12_ROOT_PARAMETER_TYPE_CBV;
		rootParameters[0].ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL;
		rootParameters[0].Descriptor.ShaderRegister = 0;
		rootParameters[0].Descriptor.RegisterSpace = 0;


		rootParameters[1].ParameterType = D3D12_ROOT_PARAMETER_TYPE_DESCRIPTOR_TABLE;
		rootParameters[1].ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL;
		rootParameters[1].DescriptorTable.NumDescriptorRanges = 1;
		rootParameters[1].DescriptorTable.pDescriptorRanges = &descriptorRange[0];

		// サンプラーの設定(s0)
		D3D12_STATIC_SAMPLER_DESC	staticSamplerDesc = {};
		staticSamplerDesc.Filter = D3D12_FILTER_MIN_MAG_MIP_LINEAR;
		staticSamplerDesc.AddressU = D3D12_TEXTURE_ADDRESS_MODE_WRAP;
		staticSamplerDesc.AddressV = D3D12_TEXTURE_ADDRESS_MODE_WRAP;
		staticSamplerDesc.AddressW = D3D12_TEXTURE_ADDRESS_MODE_WRAP;
		staticSamplerDesc.MipLODBias = 0.0f;
		staticSamplerDesc.MaxAnisotropy = 16;
		staticSamplerDesc.ComparisonFunc = D3D12_COMPARISON_FUNC_NEVER;
		staticSamplerDesc.BorderColor = D3D12_STATIC_BORDER_COLOR_TRANSPARENT_BLACK;
		staticSamplerDesc.MinLOD = 0.0f;
		staticSamplerDesc.MaxLOD = D3D12_FLOAT32_MAX;
		staticSamplerDesc.ShaderRegister = 0;
		staticSamplerDesc.RegisterSpace = 0;
		staticSamplerDesc.ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL;

		//ルートシグネチャーを作成
		D3D12_ROOT_SIGNATURE_DESC	rootSignatureDesc;
		rootSignatureDesc.NumParameters = (UINT)std::size(rootParameters);
		rootSignatureDesc.pParameters = rootParameters;
		rootSignatureDesc.NumStaticSamplers = 1;
		rootSignatureDesc.pStaticSamplers = &staticSamplerDesc;
		rootSignatureDesc.Flags = D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT;

		ComPtr<ID3DBlob>	blob = {};
		ComPtr<ID3DBlob>	error;
		if (FAILED(D3D12SerializeRootSignature(&rootSignatureDesc, D3D_ROOT_SIGNATURE_VERSION_1, &blob, &error))) {
			return FALSE;
		}
		if (FAILED(g_device->CreateRootSignature(0, blob->GetBufferPointer(), blob->GetBufferSize(), IID_PPV_ARGS(&g_rootSignature)))) {
			return FALSE;
		}
	}

	// シェーダーをコンパイル
	{
		ComPtr<ID3DBlob>	vertexShader;
		ComPtr<ID3DBlob>	pixelShader;

#if defined(_DEBUG)
		// グラフィックデバッグツールによるシェーダーのデバッグを有効にする
		UINT	compileFlags = D3DCOMPILE_DEBUG | D3DCOMPILE_SKIP_OPTIMIZATION;
#else
		UINT	compileFlags = 0;
#endif

		if (FAILED(D3DCompileFromFile(L"../Source/shaders.hlsl", nullptr, nullptr, "VSMain", "vs_5_0", compileFlags, 0, &vertexShader, nullptr))) {
			return FALSE;
		}
		if (FAILED(D3DCompileFromFile(L"../Source/shaders.hlsl", nullptr, nullptr, "PSMain", "ps_5_0", compileFlags, 0, &pixelShader, nullptr))) {
			return FALSE;
		}

		// 頂点入力レイアウトを定義
		D3D12_INPUT_ELEMENT_DESC	inputElementDescs[] = {
			{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT,		0,  0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
			{ "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT,			0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 }
		};

		// グラフィックスパイプラインの状態オブジェクトを作成
		D3D12_GRAPHICS_PIPELINE_STATE_DESC	psoDesc = {};
		psoDesc.InputLayout = { inputElementDescs, static_cast<UINT>(std::size(inputElementDescs)) };
		psoDesc.pRootSignature = g_rootSignature.Get();
		{
			D3D12_SHADER_BYTECODE	shaderBytecode;
			shaderBytecode.pShaderBytecode = vertexShader->GetBufferPointer();
			shaderBytecode.BytecodeLength = vertexShader->GetBufferSize();
			psoDesc.VS = shaderBytecode;
		}
		{
			D3D12_SHADER_BYTECODE	shaderBytecode;
			shaderBytecode.pShaderBytecode = pixelShader->GetBufferPointer();
			shaderBytecode.BytecodeLength = pixelShader->GetBufferSize();
			psoDesc.PS = shaderBytecode;
		}
		{
			D3D12_RASTERIZER_DESC	rasterizerDesc = {};
			rasterizerDesc.FillMode = D3D12_FILL_MODE_SOLID;
			rasterizerDesc.CullMode = D3D12_CULL_MODE_BACK;
			rasterizerDesc.FrontCounterClockwise = FALSE;
			rasterizerDesc.DepthBias = D3D12_DEFAULT_DEPTH_BIAS;
			rasterizerDesc.DepthBiasClamp = D3D12_DEFAULT_DEPTH_BIAS_CLAMP;
			rasterizerDesc.SlopeScaledDepthBias = D3D12_DEFAULT_SLOPE_SCALED_DEPTH_BIAS;
			rasterizerDesc.DepthClipEnable = TRUE;
			rasterizerDesc.MultisampleEnable = FALSE;
			rasterizerDesc.AntialiasedLineEnable = FALSE;
			rasterizerDesc.ForcedSampleCount = 0;
			rasterizerDesc.ConservativeRaster = D3D12_CONSERVATIVE_RASTERIZATION_MODE_OFF;
			psoDesc.RasterizerState = rasterizerDesc;
		}
		{
			D3D12_BLEND_DESC	blendDesc = {};
			blendDesc.AlphaToCoverageEnable = FALSE;
			blendDesc.IndependentBlendEnable = FALSE;
			for (UINT i = 0; i < D3D12_SIMULTANEOUS_RENDER_TARGET_COUNT; ++i) {
				blendDesc.RenderTarget[i].BlendEnable = FALSE;
				blendDesc.RenderTarget[i].LogicOpEnable = FALSE;
				blendDesc.RenderTarget[i].SrcBlend = D3D12_BLEND_ONE;
				blendDesc.RenderTarget[i].DestBlend = D3D12_BLEND_ZERO;
				blendDesc.RenderTarget[i].BlendOp = D3D12_BLEND_OP_ADD;
				blendDesc.RenderTarget[i].SrcBlendAlpha = D3D12_BLEND_ONE;
				blendDesc.RenderTarget[i].DestBlendAlpha = D3D12_BLEND_ZERO;
				blendDesc.RenderTarget[i].BlendOpAlpha = D3D12_BLEND_OP_ADD;
				blendDesc.RenderTarget[i].LogicOp = D3D12_LOGIC_OP_NOOP;
				blendDesc.RenderTarget[i].RenderTargetWriteMask = D3D12_COLOR_WRITE_ENABLE_ALL;
			}
			psoDesc.BlendState = blendDesc;
		}
		psoDesc.DepthStencilState.DepthEnable = FALSE;
		psoDesc.DepthStencilState.StencilEnable = FALSE;
		psoDesc.SampleMask = UINT_MAX;
		psoDesc.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
		psoDesc.NumRenderTargets = 1;
		psoDesc.RTVFormats[0] = DXGI_FORMAT_R8G8B8A8_UNORM;
		psoDesc.SampleDesc.Count = 1;

		if (FAILED(g_device->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&g_pipelineState)))) {
			return FALSE;
		}
	}

	// コマンドリストを作成
	if (FAILED(g_device->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, g_commandAllocator.Get(), g_pipelineState.Get(), IID_PPV_ARGS(&g_commandList)))) {
		return FALSE;
	}

	if (FAILED(g_commandList->Close())) {
		return FALSE;
	}

	// 頂点バッファを作成
	{
		//#TODO 三角形一つで画面いっぱいに描く、というのもやってみたい　uvなどがめんどくさそうだけど
		// ジオメトリを定義
		Vertex	triangleVertices[] = {
			{ { -1.0f,  1.0f * g_aspectRatio, 0.0f },{ 0.0f, 0.0f } },
			{ {  1.0f,  1.0f * g_aspectRatio, 0.0f },{ 1.0f, 0.0f } },
			{ {  1.0f, -1.0f * g_aspectRatio, 0.0f },{ 1.0f, 1.0f } },
			{ { -1.0f,  1.0f * g_aspectRatio, 0.0f },{ 0.0f, 0.0f } },
			{ {  1.0f, -1.0f * g_aspectRatio, 0.0f },{ 1.0f, 1.0f } },
			{ { -1.0f, -1.0f * g_aspectRatio, 0.0f },{ 0.0f, 1.0f } },
		};

		const UINT	vertexBufferSize = sizeof(triangleVertices);

		{
			D3D12_HEAP_PROPERTIES	heapProperties = {};
			heapProperties.Type = D3D12_HEAP_TYPE_UPLOAD;
			heapProperties.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
			heapProperties.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
			heapProperties.CreationNodeMask = 1;
			heapProperties.VisibleNodeMask = 1;

			D3D12_RESOURCE_DESC	resourceDesc = {};
			resourceDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
			resourceDesc.Alignment = 0;
			resourceDesc.Width = vertexBufferSize;
			resourceDesc.Height = 1;
			resourceDesc.DepthOrArraySize = 1;
			resourceDesc.MipLevels = 1;
			resourceDesc.Format = DXGI_FORMAT_UNKNOWN;
			resourceDesc.SampleDesc.Count = 1;
			resourceDesc.SampleDesc.Quality = 0;
			resourceDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
			resourceDesc.Flags = D3D12_RESOURCE_FLAG_NONE;

			if (FAILED(g_device->CreateCommittedResource(&heapProperties, D3D12_HEAP_FLAG_NONE, &resourceDesc, D3D12_RESOURCE_STATE_GENERIC_READ, nullptr, IID_PPV_ARGS(&g_vertexBuffer)))) {
				return FALSE;
			}
		}

		// 頂点バッファに頂点データをコピー
		UINT8*		pVertexDataBegin = nullptr;
		D3D12_RANGE	readRange = { 0, 0 };		// CPUからバッファを読み込まない設定 nullptrでもよさそう？
		if (FAILED(g_vertexBuffer->Map(0, &readRange, reinterpret_cast<void**>(&pVertexDataBegin)))) {
			return FALSE;
		}
		memcpy(pVertexDataBegin, triangleVertices, sizeof(triangleVertices));
		g_vertexBuffer->Unmap(0, nullptr);

		// 頂点バッファビューを初期化
		g_vertexBufferView.BufferLocation = g_vertexBuffer->GetGPUVirtualAddress();
		g_vertexBufferView.StrideInBytes = sizeof(Vertex);
		g_vertexBufferView.SizeInBytes = vertexBufferSize;

		////平面のインデックス
		//const UINT	vertexBufferSize = sizeof(triangleVertices);


		////頂点バッファの作成
		//resource_desc.Width = sizeof(Vertex3D) * VERT_NUM * ARC_NUM;
		//hr = device->CreateCommittedResource(&heap_properties, D3D12_HEAP_FLAG_NONE, &resource_desc, D3D12_RESOURCE_STATE_GENERIC_READ, nullptr, IID_PPV_ARGS(&vertex_buffer_));
		//if (FAILED(hr)) {
		//	return hr;
		//}


		////平面のバッファ
		//UINT16* pIndexDataBegin = nullptr;
		//if (FAILED(g_indexBuffer->Map(0, &readRange, reinterpret_cast<void**>(&pIndexDataBegin)))) {
		//	return FALSE;
		//}
		//memcpy(pIndexDataBegin, , sizeof());
		//g_indexBuffer->Unmap(0, nullptr);

		////インデックスバッファ
		//g_indexBufferView.BufferLocation = g_indexBuffer->GetGPUVirtualAddress();
		//g_indexBufferView.SizeInBytes = sizeof(planeIndices);
		//g_indexBufferView.Format = DXGI_FORMAT_R16_UINT;
	}

	//定数バッファの初期化
	{
		D3D12_HEAP_PROPERTIES	heapProperties = {};
		heapProperties.Type = D3D12_HEAP_TYPE_UPLOAD;
		heapProperties.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
		heapProperties.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
		heapProperties.CreationNodeMask = 1;
		heapProperties.VisibleNodeMask = 1;

		D3D12_RESOURCE_DESC	resourceDesc = {};
		resourceDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
		resourceDesc.Alignment = 0;
		resourceDesc.Width = 256;
		resourceDesc.Height = 1;
		resourceDesc.DepthOrArraySize = 1;
		resourceDesc.MipLevels = 1;
		resourceDesc.Format = DXGI_FORMAT_UNKNOWN;
		resourceDesc.SampleDesc.Count = 1;
		resourceDesc.SampleDesc.Quality = 0;
		resourceDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
		resourceDesc.Flags = D3D12_RESOURCE_FLAG_NONE;

		if (FAILED(g_device->CreateCommittedResource(&heapProperties, D3D12_HEAP_FLAG_NONE, &resourceDesc, D3D12_RESOURCE_STATE_GENERIC_READ, nullptr, IID_PPV_ARGS(&g_constantBuffer)))) {
			return FALSE;
		}
	}

	// テクスチャ関係の初期化
	if (!InitTexture(hWnd)) {
		return FALSE;
	}

	return TRUE;
}

// テクスチャを初期化
BOOL InitTexture(HWND hWnd)
{
	{
		// ビットマップファイルを読み込み
		// (DDB＆GetPixelでピクセルデータを読み込む簡易版なのでアルファ値には対応してません)
		HBITMAP	hBitmap;
		hBitmap = (HBITMAP)LoadImage(0, _T("../Resource/texture.bmp"), IMAGE_BITMAP, 0, 0, LR_LOADFROMFILE);
		if (!hBitmap) {
			return FALSE;
		}

		HDC		hMemDC;
		BITMAP	bmp;
		hMemDC = CreateCompatibleDC(NULL);
		SelectObject(hMemDC, hBitmap);
		GetObject(hBitmap, sizeof(BITMAP), &bmp);
		DeleteObject(hBitmap);
		std::vector<uint32_t>	pixel;
		for (int y = 0; y < bmp.bmHeight; y++) {
			for (int x = 0; x < bmp.bmWidth; x++) {
				COLORREF	color = GetPixel(hMemDC, x, y);
				BYTE		r = color & 0xFF;
				BYTE		g = (color >> 8) & 0xFF;
				BYTE		b = (color >> 16) & 0xFF;
				pixel.push_back((r << 16) | (g << 8) | b);
			}
		}
		DeleteObject(hMemDC);

		//テクスチャ用のリソースを作成
		D3D12_HEAP_PROPERTIES	heapProperties = {};
		heapProperties.Type = D3D12_HEAP_TYPE_CUSTOM;
		heapProperties.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_WRITE_BACK;
		heapProperties.MemoryPoolPreference = D3D12_MEMORY_POOL_L0;
		heapProperties.CreationNodeMask = 1;
		heapProperties.VisibleNodeMask = 1;

		D3D12_RESOURCE_DESC		resourceDesc = {};
		resourceDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
		resourceDesc.Width = bmp.bmWidth;
		resourceDesc.Height = bmp.bmHeight;
		resourceDesc.DepthOrArraySize = 1;
		resourceDesc.MipLevels = 1;
		resourceDesc.Format = DXGI_FORMAT_B8G8R8A8_UNORM;
		resourceDesc.SampleDesc.Count = 1;
		resourceDesc.SampleDesc.Quality = 0;
		resourceDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
		resourceDesc.Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN;
		if (FAILED(g_device->CreateCommittedResource(&heapProperties, D3D12_HEAP_FLAG_NONE, &resourceDesc, D3D12_RESOURCE_STATE_GENERIC_READ, nullptr, IID_PPV_ARGS(&g_rTexture)))) {
			return FALSE;
		}

		// テクスチャー用の記述子ヒープを作成
		D3D12_DESCRIPTOR_HEAP_DESC	descriptorHeapDesc = {};
		descriptorHeapDesc.NumDescriptors = 1 + 1;	//t0,t1ということ　わかりづらいのでいい感じにしたい
		descriptorHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
		descriptorHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
		descriptorHeapDesc.NodeMask = 0;
		if (FAILED(g_device->CreateDescriptorHeap(&descriptorHeapDesc, IID_PPV_ARGS(&g_dhTexture)))) {
			return FALSE;
		}

		// テクスチャー用のシェーダーリソースビューを作成
		D3D12_CPU_DESCRIPTOR_HANDLE		cpuDescriptorHandle = {};
		D3D12_SHADER_RESOURCE_VIEW_DESC	shaderResourceViewDesc = {};
		shaderResourceViewDesc.Format = DXGI_FORMAT_B8G8R8A8_UNORM;
		shaderResourceViewDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
		shaderResourceViewDesc.Texture2D.MipLevels = 1;
		shaderResourceViewDesc.Texture2D.MostDetailedMip = 0;
		shaderResourceViewDesc.Texture2D.PlaneSlice = 0;
		shaderResourceViewDesc.Texture2D.ResourceMinLODClamp = 0.0F;
		shaderResourceViewDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
		cpuDescriptorHandle = g_dhTexture->GetCPUDescriptorHandleForHeapStart();
		g_device->CreateShaderResourceView(g_rTexture.Get(), &shaderResourceViewDesc, cpuDescriptorHandle);

		// 画像データをサブリソースへコピー
		D3D12_BOX	box = { 0, 0, 0, (UINT)bmp.bmWidth, (UINT)bmp.bmHeight, 1 };
		if (FAILED(g_rTexture->WriteToSubresource(0, &box, &pixel[0], sizeof(uint32_t)*bmp.bmWidth, sizeof(uint32_t)*bmp.bmWidth*bmp.bmHeight))) {
			return FALSE;
		}
	}

	// 同期オブジェクトを作成してリソースがGPUにアップロードされるまで待機
	{
		if (FAILED(g_device->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&g_fence)))) return FALSE;
		g_fenceValue = 1;

		// フレーム同期に使用するイベントハンドラを作成
		g_fenceEvent = CreateEvent(nullptr, FALSE, FALSE, nullptr);
		if (g_fenceEvent == nullptr) {
			if (FAILED(HRESULT_FROM_WIN32(GetLastError()))) return FALSE;
		}

		if (!WaitForPreviousFrame()) return FALSE;
	}

	//環境マップ読み込み
	{
		wchar_t ddsFileName[] = L"../Resource/uffizi.dds";
		TexMetadata mdata;
		HRESULT hr = GetMetadataFromDDSFile(ddsFileName, DDS_FLAGS_NONE, mdata);
		if (FAILED(hr))
		{
			//swprintf_s(buff, L"ファイルの読み込みに失敗しました　Failed to open texture file\n\nFilename = %ls\nHRESULT %08X", lpCmdLine, hr);
			MessageBox(hWnd, _T("ファイルの読み込みに失敗しました"), _T("GetMetadataFromDDSFile"), MB_OK | MB_ICONEXCLAMATION);
			return FALSE;
		}

		switch (mdata.format)
		{
		case DXGI_FORMAT_BC6H_TYPELESS:
		case DXGI_FORMAT_BC6H_UF16:
		case DXGI_FORMAT_BC6H_SF16:
		case DXGI_FORMAT_BC7_TYPELESS:
		case DXGI_FORMAT_BC7_UNORM:
		case DXGI_FORMAT_BC7_UNORM_SRGB:
			if (g_featureLevel < D3D_FEATURE_LEVEL_11_0)
			{
				wchar_t buff[2048] = {};
				swprintf_s(buff, L"BC6H/BC7 requires DirectX 11 hardware\n\nFilename = %ls\nDXGI Format %d\nFeature Level %d", ddsFileName, mdata.format, g_featureLevel);
				MessageBoxW(nullptr, buff, L"DDSView", MB_OK | MB_ICONEXCLAMATION);
				return 0;
			}
			break;

		default:
		{
			//DX12にはないみたい
			//UINT flags = 0;
			//hr = g_device->CheckFormatSupport(mdata.format, &flags);
			//if (FAILED(hr) || !(flags & (D3D11_FORMAT_SUPPORT_TEXTURE1D | D3D11_FORMAT_SUPPORT_TEXTURE2D | D3D11_FORMAT_SUPPORT_TEXTURE3D)))
			//{
			//	wchar_t buff[2048] = {};
			//	swprintf_s(buff, L"Format not supported by DirectX hardware\n\nFilename = %ls\nDXGI Format %d\nFeature Level %d\nHRESULT = %08X", lpCmdLine, mdata.format, g_featureLevel, hr);
			//	MessageBoxW(nullptr, buff, L"DDSView", MB_OK | MB_ICONEXCLAMATION);
			//	return 0;
			//}
		}
		break;
		}

		ScratchImage image;
		hr = LoadFromDDSFile(ddsFileName, DDS_FLAGS_NONE, &mdata, image);
		if (FAILED(hr))
		{
			wchar_t buff[2048] = {};
			swprintf_s(buff, L"Failed to load texture file\n\nFilename = %ls\nHRESULT %08X", ddsFileName, hr);
			MessageBoxW(nullptr, buff, L"DDSView", MB_OK | MB_ICONEXCLAMATION);
			return 0;
		}

		//テクスチャ用のリソースを作成
		D3D12_HEAP_PROPERTIES	heapProperties = {};
		heapProperties.Type = D3D12_HEAP_TYPE_CUSTOM;
		heapProperties.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_WRITE_BACK;
		heapProperties.MemoryPoolPreference = D3D12_MEMORY_POOL_L0;
		heapProperties.CreationNodeMask = 1;
		heapProperties.VisibleNodeMask = 1;

		D3D12_RESOURCE_DESC		resourceDesc = {};
		resourceDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
		resourceDesc.Width = mdata.width;
		resourceDesc.Height = (UINT)mdata.height;
		resourceDesc.DepthOrArraySize = 6;
		resourceDesc.MipLevels = 1;
		resourceDesc.Format = DXGI_FORMAT_R16G16B16A16_FLOAT;
		resourceDesc.SampleDesc.Count = 1;
		resourceDesc.SampleDesc.Quality = 0;
		resourceDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
		resourceDesc.Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN;
		if (FAILED(g_device->CreateCommittedResource(&heapProperties, D3D12_HEAP_FLAG_NONE, &resourceDesc, D3D12_RESOURCE_STATE_GENERIC_READ, nullptr, IID_PPV_ARGS(&g_rTestTexture)))) {
			return FALSE;
		}

		// テクスチャー用のシェーダーリソースビューを作成
		//#TODO そのまま書いているので、もう少し良い感じに
		D3D12_CPU_DESCRIPTOR_HANDLE		cpuDescriptorHandle = {};
		D3D12_SHADER_RESOURCE_VIEW_DESC	shaderResourceViewDesc = {};
		shaderResourceViewDesc.Format = DXGI_FORMAT_R16G16B16A16_FLOAT;
		shaderResourceViewDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURECUBE;
		//shaderResourceViewDesc.TextureCubeArray.NumCubes = resourceDesc.DepthOrArraySize / 6;
		shaderResourceViewDesc.Texture2D.MipLevels = 1;
		shaderResourceViewDesc.Texture2D.MostDetailedMip = 0;
		shaderResourceViewDesc.Texture2D.PlaneSlice = 0;
		shaderResourceViewDesc.Texture2D.ResourceMinLODClamp = 0.0F;
		shaderResourceViewDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
		cpuDescriptorHandle = g_dhTexture->GetCPUDescriptorHandleForHeapStart();
		cpuDescriptorHandle.ptr += g_device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);;
		g_device->CreateShaderResourceView(g_rTestTexture.Get(), &shaderResourceViewDesc, cpuDescriptorHandle);

		// 画像データをサブリソースへコピー
		//UpdateSubresourcesでやったら簡単？
		for (int i = 0; i < 6; i++) {
			const Image* pImage = image.GetImage(0, i, 0);
			D3D12_BOX	box = { 0, 0, 0, (UINT)pImage->width, (UINT)pImage->height, 1 };
			if (FAILED(g_rTestTexture->WriteToSubresource(i, &box, pImage->pixels, (UINT)pImage->rowPitch, (UINT)pImage->slicePitch))) {
				return FALSE;
			}
		}

		//g_texCommandList->Close();
		//ID3D12CommandList	*ppCommandLists[] = { g_texCommandList.Get() };
		//g_commandQueue->ExecuteCommandLists(_countof(ppCommandLists), ppCommandLists);
	}

	// 同期オブジェクトを作成してリソースがGPUにアップロードされるまで待機
	{
		if (FAILED(g_device->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&g_fence)))) return FALSE;
		g_fenceValue = 1;

		// フレーム同期に使用するイベントハンドラを作成
		g_fenceEvent = CreateEvent(nullptr, FALSE, FALSE, nullptr);
		if (g_fenceEvent == nullptr) {
			if (FAILED(HRESULT_FROM_WIN32(GetLastError()))) return FALSE;
		}

		if (!WaitForPreviousFrame()) return FALSE;
	}














	return TRUE;
}

// 描画
BOOL Draw()
{
	if (FAILED(g_commandAllocator->Reset())) return FALSE;
	if (FAILED(g_commandList->Reset(g_commandAllocator.Get(), g_pipelineState.Get()))) return FALSE;

	g_commandList->SetGraphicsRootSignature(g_rootSignature.Get());
	g_commandList->RSSetViewports(1, &g_viewport);
	g_commandList->RSSetScissorRects(1, &g_scissorRect);
	
	//定数バッファをシェーダーのレジスタへ設定
	{
		float* buffer = {};
		HRESULT hr = g_constantBuffer->Map(0, nullptr, (void**)&buffer);
		if (FAILED(hr)) {
			return hr;
		}

		buffer[0] = g_time;

		g_constantBuffer->Unmap(0, nullptr);
		buffer = nullptr;

		g_commandList->SetGraphicsRootConstantBufferView(0, g_constantBuffer->GetGPUVirtualAddress());

		g_time += 0.01f;
	}

	// テクスチャをシェーダーのレジスタへ設定
	{
		ID3D12DescriptorHeap* descriptorHeap[] = { g_dhTexture.Get() };
		g_commandList->SetDescriptorHeaps((UINT)std::size(descriptorHeap), descriptorHeap);
		g_commandList->SetGraphicsRootDescriptorTable(1, g_dhTexture->GetGPUDescriptorHandleForHeapStart());
	}

	// バックバッファをレンダリングターゲットとして使用
	{
		D3D12_RESOURCE_BARRIER	resourceBarrier = {};
		resourceBarrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
		resourceBarrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
		resourceBarrier.Transition.pResource = g_renderTargets[g_frameIndex].Get();
		resourceBarrier.Transition.StateBefore = D3D12_RESOURCE_STATE_PRESENT;
		resourceBarrier.Transition.StateAfter = D3D12_RESOURCE_STATE_RENDER_TARGET;
		resourceBarrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
		g_commandList->ResourceBarrier(1, &resourceBarrier);
	}

	D3D12_CPU_DESCRIPTOR_HANDLE	rtvHandle = {};
	rtvHandle.ptr = g_descriptorHeap->GetCPUDescriptorHandleForHeapStart().ptr + g_frameIndex * g_descriptorSize;
	g_commandList->OMSetRenderTargets(1, &rtvHandle, FALSE, nullptr);

	// バックバッファに描画
	const float	clearColor[] = { 0.0f, 0.2f, 0.4f, 1.0f };
	g_commandList->ClearRenderTargetView(rtvHandle, clearColor, 0, nullptr);
	g_commandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
	g_commandList->IASetVertexBuffers(0, 1, &g_vertexBufferView);
	//g_commandList->IASetIndexBuffer(&g_indexView);
	g_commandList->DrawInstanced(6, 2, 0, 0);

	// バックバッファを表示
	{
		D3D12_RESOURCE_BARRIER	resourceBarrier = {};
		resourceBarrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
		resourceBarrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
		resourceBarrier.Transition.pResource = g_renderTargets[g_frameIndex].Get();
		resourceBarrier.Transition.StateBefore = D3D12_RESOURCE_STATE_RENDER_TARGET;
		resourceBarrier.Transition.StateAfter = D3D12_RESOURCE_STATE_PRESENT;
		resourceBarrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
		g_commandList->ResourceBarrier(1, &resourceBarrier);
	}

	if (FAILED(g_commandList->Close())) return FALSE;


	// コマンドリストを実行
	ID3D12CommandList	*ppCommandLists[] = { g_commandList.Get() };
	g_commandQueue->ExecuteCommandLists((UINT)std::size(ppCommandLists), ppCommandLists);

	// フレームを最終出力
	if (FAILED(g_swapChain->Present(1, 0))) return FALSE;

	return WaitForPreviousFrame();
}

BOOL WaitForPreviousFrame()
{
	const UINT64	fence = g_fenceValue;
	if (FAILED(g_commandQueue->Signal(g_fence.Get(), fence))) return FALSE;
	g_fenceValue++;

	// 前のフレームが終了するまで待機
	if (g_fence->GetCompletedValue() < fence) {
		if (FAILED(g_fence->SetEventOnCompletion(fence, g_fenceEvent))) return FALSE;
		WaitForSingleObject(g_fenceEvent, INFINITE);
	}

	g_frameIndex = g_swapChain->GetCurrentBackBufferIndex();

	return TRUE;
}
