#include "APP.h"

App *g_App = 0;
//------------------------------------------------------------------------------
//static callback function
//------------------------------------------------------------------------------
static void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	g_App->onKey(window, key, scancode, action, mods);
}

static void mousebutton_callback(GLFWwindow* window, int button, int action, int mods)
{
	g_App->onMouseButton(window, button, action, mods);
}

static void cursorpos_callback(GLFWwindow* window, double x, double y)
{
	g_App->onMouseMove(window, x, y);
}

static void resize_callback(GLFWwindow* window, int w, int h)
{
	g_App->onResize(window, w, h);
}

static void mousewheel_callback(GLFWwindow* window, double x, double y)
{
	g_App->onMouseWheel(window, x, y);
}
////////////////////////////////////////////////////////////////////////////////////
App::App() :mTitle("sample app"), mWidth(1024),
mHeight(768), isTransparency(false), isLineSmooth(false)
{
	g_App = this;
}
App::App(int w, int h) : mTitle("sample app"), mWidth(w),
mHeight(h), isTransparency(false), isLineSmooth(false)
{
	g_App = this;
}
App::~App()
{
}

float App::AspectRatio()const
{
	return static_cast<float>(mWidth) / mHeight;
}

int App::Run()
{
	// GLFW rendering process
	while (!glfwWindowShouldClose(window))
	{
		// Display routine
		UpdateScene();
		Rendering();
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwDestroyWindow(window);
	glfwTerminate();

	return EXIT_SUCCESS;
}

bool App::Init()
{
	if (!InitGLFW())
		return false;

	// callbacks initialization
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mousebutton_callback);
	glfwSetCursorPosCallback(window, cursorpos_callback);
	glfwSetScrollCallback(window, mousewheel_callback);
	glfwSetWindowSizeCallback(window, resize_callback);

	// Initialize GLEW
	glewExperimental = true; // Needed in core profile
	if (glewInit() != GLEW_OK) {
		cout << "Failed to initialize GLEW" << endl;
		return false;
	}

	if (!InitGL())
		return false;

	return true;
}

bool App::InitGLFW()
{
	// GLFW initialization
	glfwSetErrorCallback(error_callback);
	if (!glfwInit()) return EXIT_FAILURE;

	window = glfwCreateWindow(mWidth, mHeight, mTitle.c_str(), NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return EXIT_FAILURE;
	}

	return true;
}

bool App::InitGL()
{
	/*
	// polygon fill mode
	glPolygonMode(GL_BACK, GL_FILL);

	// enable depth testing
	//glEnable(GL_DEPTH_TEST);
	//glDepthFunc(GL_LESS);

	glEnable(GL_NORMALIZE);

	// transparency settings
	if (isTransparency)
	{
	glEnable(GL_ALPHA_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}

	glEnable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_POLYGON_OFFSET_LINE);
	glPolygonOffset((float) 1.0, (float) 1e-5);

	// line anti-aliasing
	if (isLineSmooth)
	{
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	}

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	*/
	return true;
}
void App::onResize(GLFWwindow* window, int w, int h)
{
	mWidth = w;
	mHeight = h;
}
void App::onKey(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// ESC
	if ((key == GLFW_KEY_ESCAPE) && (action == GLFW_PRESS))
	{
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
}
void App::onMouseWheel(GLFWwindow* window, double x, double y)
{
}
void App::onMouseMove(GLFWwindow* window, double xd, double yd)
{
}
void App::onMouseButton(GLFWwindow* window, int button, int action, int mods)
{
}
void App::UpdateScene()
{
}
void App::Rendering()
{
}