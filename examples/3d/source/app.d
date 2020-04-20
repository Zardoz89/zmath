import std.stdio;
import std.string : fromStringz;
import std.conv;

import bindbc.sdl;
import bindbc.opengl;

import zmath;

enum W_WIDTH = 800;
enum W_HEIGHT = 600;

SDL_Window* appWin;
SDL_GLContext glContext;
SDL_Surface* imgSurf;
GLuint textureID = 0;
GLuint vertexArrayId;
GLuint vertexBufferId;
GLuint colorBufferId;

const string VertexShaderCode = "
#version 330 core

// Input vertex data
layout(location = 0) in vec3 vertexPosition_modelspace;
// Notice that the '1' here equals the '1' in glVertexAttribPointer
layout(location = 1) in vec3 vertexColor;

// ModelViewProyection matrix
uniform mat4 MVP;

out vec3 fragmentColor;

void main(){
  gl_Position = MVP * vec4(vertexPosition_modelspace, 1);
  fragmentColor = vertexColor;
}
";
const string FragmentShaderCode = "
#version 330 core

// Interpolated values from the vertex shaders
in vec3 fragmentColor;

out vec3 color;

void main(){
  color = fragmentColor;
}
";

static const GLfloat[12*3*3] g_vertex_buffer_data = [
    -1.0f,-1.0f,-1.0f, // triangle 1 : begin
    -1.0f,-1.0f, 1.0f,
    -1.0f, 1.0f, 1.0f, // triangle 1 : end
    1.0f, 1.0f,-1.0f, // triangle 2 : begin
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f,-1.0f, // triangle 2 : end
    1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f,-1.0f,
    1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    -1.0f,-1.0f, 1.0f,
    1.0f,-1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f,-1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f, 1.0f,-1.0f,
    -1.0f, 1.0f,-1.0f,
    1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f, 1.0f,
    1.0f,-1.0f, 1.0f
];

static const GLfloat[12*3*3] g_color_buffer_data = [
    0.583f,  0.771f,  0.014f,
    0.609f,  0.115f,  0.436f,
    0.327f,  0.483f,  0.844f,
    0.822f,  0.569f,  0.201f,
    0.435f,  0.602f,  0.223f,
    0.310f,  0.747f,  0.185f,
    0.597f,  0.770f,  0.761f,
    0.559f,  0.436f,  0.730f,
    0.359f,  0.583f,  0.152f,
    0.483f,  0.596f,  0.789f,
    0.559f,  0.861f,  0.639f,
    0.195f,  0.548f,  0.859f,
    0.014f,  0.184f,  0.576f,
    0.771f,  0.328f,  0.970f,
    0.406f,  0.615f,  0.116f,
    0.676f,  0.977f,  0.133f,
    0.971f,  0.572f,  0.833f,
    0.140f,  0.616f,  0.489f,
    0.997f,  0.513f,  0.064f,
    0.945f,  0.719f,  0.592f,
    0.543f,  0.021f,  0.978f,
    0.279f,  0.317f,  0.505f,
    0.167f,  0.620f,  0.077f,
    0.347f,  0.857f,  0.137f,
    0.055f,  0.953f,  0.042f,
    0.714f,  0.505f,  0.345f,
    0.783f,  0.290f,  0.734f,
    0.722f,  0.645f,  0.174f,
    0.302f,  0.455f,  0.848f,
    0.225f,  0.587f,  0.040f,
    0.517f,  0.713f,  0.338f,
    0.053f,  0.959f,  0.120f,
    0.393f,  0.621f,  0.362f,
    0.673f,  0.211f,  0.457f,
    0.820f,  0.883f,  0.371f,
    0.982f,  0.099f,  0.879f
];


/// Initialize SDL2
bool initSDL()
{
  // Load SDL libs
  const SDLSupport ret = loadSDL();
  if(ret != sdlSupport) {
    writeln("Error loading SDL dll");
    return false;
  }

  // Initialise SDL
  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    writeln("SDL_Init: ", fromStringz(SDL_GetError()));
    return false;
  }
  return true;
}

/// Set OpenGL flags to start a OpenGL 3 context
void setGlFlags()
{
  version(OSX) {
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG); // Always required on Mac
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
  } else {
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
  }
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
}

/// Creates the SDL2 window
bool createAppWindow()
{
  // Create a window
  enum windowFlags = SDL_WINDOW_OPENGL | SDL_WINDOW_ALLOW_HIGHDPI | SDL_WINDOW_SHOWN;
  appWin = SDL_CreateWindow(
      "Camera example",
      SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
      W_WIDTH, W_HEIGHT,
      windowFlags
      );
  if (appWin is null) {
    writefln("SDL_CreateWindow: ", fromStringz(SDL_GetError()));
    return false;
  }
  return true;
}

/// Initialize OpenGL stuff
bool initGl()
{
  glContext = SDL_GL_CreateContext(appWin);
  if (glContext == null) {
    writeln("OpenGL context couldn't be created! SDL Error: ", fromStringz(SDL_GetError()));
    return false;
  }

  const GLSupport openglLoaded = loadOpenGL();
  if ( openglLoaded != glSupport) {
    writeln("Error loading OpenGL shared library", to!string(openglLoaded));
    return false;
  }
  SDL_GL_MakeCurrent(appWin, glContext);

  SDL_GL_SetSwapInterval(1); // Enable VSync
  glEnable(GL_TEXTURE_2D); // Enable textures
  glEnable(GL_DEPTH_TEST); // Enable depth test
  glDepthFunc(GL_LESS); // Accept fragment if it closer to the camera than the former one

  // Clear buffers
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClearColor( 0f, 0f, 0f, 1f );

  return true;
}

GLuint loadShaders()
{
  import std.conv : to;
  import std.string : toStringz;

  // Show compiling or linking error
  void showError(GLuint id, int infoLogLength) {
    import std.string : fromStringz;
    char[] errorMessage;
    errorMessage.length = infoLogLength + 1;
    glGetShaderInfoLog(id, infoLogLength, null, errorMessage.ptr);
    writeln(fromStringz(errorMessage.ptr));
  }

  GLuint vertexShaderId = glCreateShader(GL_VERTEX_SHADER);
  GLuint fragmentShaderId = glCreateShader(GL_FRAGMENT_SHADER);
  GLint result = GL_FALSE;
  int infoLogLength;

  // Compile Vertex Shader
  writeln("Compiling vertex shader");
  const char* vertexShaderCodePtr = toStringz(VertexShaderCode);
  glShaderSource(vertexShaderId, 1, &vertexShaderCodePtr, null);
  glCompileShader(vertexShaderId);

  // check Vertex Shader
  glGetShaderiv(vertexShaderId, GL_COMPILE_STATUS, &result);
  glGetShaderiv(vertexShaderId, GL_INFO_LOG_LENGTH, &infoLogLength);
  if ( infoLogLength > 0 ){
    showError(vertexShaderId, infoLogLength);
  }

  // Compile Fragment Shader
  writeln("Compiling fragment shader");
  const char* fragmentShaderCodePtr = toStringz(FragmentShaderCode);
  glShaderSource(fragmentShaderId, 1, &fragmentShaderCodePtr, null);
  glCompileShader(fragmentShaderId);

  // check Vertex Shader
  glGetShaderiv(fragmentShaderId, GL_COMPILE_STATUS, &result);
  glGetShaderiv(fragmentShaderId, GL_INFO_LOG_LENGTH, &infoLogLength);
  if ( infoLogLength > 0 ){
    showError(fragmentShaderId, infoLogLength);
  }

  // Linking
  writeln("Linking shaders");
  GLuint programId = glCreateProgram();
  glAttachShader(programId, vertexShaderId);
  glAttachShader(programId, fragmentShaderId);
  glLinkProgram(programId);

  glGetProgramiv(programId, GL_LINK_STATUS, &result);
  glGetProgramiv(programId, GL_INFO_LOG_LENGTH, &infoLogLength);
  if ( infoLogLength > 0 ){
    showError(programId, infoLogLength);
  }

  glDetachShader(programId, vertexShaderId);
  glDetachShader(programId, fragmentShaderId);

  glDeleteShader(vertexShaderId);
  glDeleteShader(fragmentShaderId);

  return programId;
}

/// Creates a VAO and load VBOs with data of the model
void loadModel()
{
  // Generates a VAO
  glGenVertexArrays(1, &vertexArrayId);
  glBindVertexArray(vertexArrayId);

  // Fill the Vertex Buffer with vertex coord
  glGenBuffers(1, &vertexBufferId);
  glBindBuffer(GL_ARRAY_BUFFER, vertexBufferId);
  glBufferData(GL_ARRAY_BUFFER, g_vertex_buffer_data.sizeof, g_vertex_buffer_data.ptr, GL_STATIC_DRAW);

  // Fill the Vertex Buffer with color data
  glGenBuffers(1, &colorBufferId);
  glBindBuffer(GL_ARRAY_BUFFER, colorBufferId);
  glBufferData(GL_ARRAY_BUFFER, g_color_buffer_data.sizeof, g_color_buffer_data.ptr, GL_STATIC_DRAW);
}

/// Render the scene
void render()
{
  // Cleaning buffers
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClearColor( 0f, 0f, 0f, 1f );

  // Render something

  // 1st attribute buffer : vertices
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, vertexBufferId);
  glVertexAttribPointer(
      0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
      3,                  // size
      GL_FLOAT,           // type
      GL_FALSE,           // normalized?
      0,                  // stride
      null                // array buffer offset
      );

  // 2nd attribute buffer : colors
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, colorBufferId);
  glVertexAttribPointer(
      1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
      3,                                // size
      GL_FLOAT,                         // type
      GL_FALSE,                         // normalized?
      0,                                // stride
      null                              // array buffer offset
  );

  // Draw the VAO !
  glDrawArrays(GL_TRIANGLES, 0, 12*3);
  glDisableVertexAttribArray(0);

}

/// Update screen
void updateScreen()
{
  glFlush();
  SDL_GL_SwapWindow(appWin);
}

/// Closes/destroy everything
void close()
{
  if (glContext !is null) {
    SDL_GL_DeleteContext(glContext);
  }
  if (appWin !is null) {
    SDL_DestroyWindow(appWin);
  }
  SDL_Quit();
}

/// Recaulcates the MVP matrix
void recalMatrixes(bool ortho, in Qua_f cameraPos, const float modelAngle, in Vec3f modelPos, GLuint programId)
{
  const proyMatrix = proyectionMatrix(ortho);
  auto viewMatrix = translateMat(0, 0, -10f) * cast(Mat4f)cameraPos;
  writeln("view matrix:\n", viewMatrix.toString());
  auto modelMatrix = rotMat(Vec3f.Y_AXIS, modelAngle) * translateMat(modelPos);
  writeln("model matrix:\n", modelMatrix.toString());
  auto mvp = proyMatrix * viewMatrix * modelMatrix;
  writeln("mvp matrix:\n", mvp.toString());

  updateMVPMatrix(mvp, programId);
}

/// Generates a proyection matrix
auto proyectionMatrix(bool ortho)
{
  float viewPortRatio = (1.0f * W_WIDTH) / W_HEIGHT;
  Mat4f proyMat;
  if (ortho) {
    proyMat = orthoMat(10 * viewPortRatio, 10, 100f);
  } else {
    import std.math : PI_2;
    proyMat = perspectiveMat(PI_2, viewPortRatio, 0.1, 100f);
  }
  writeln("Proyection matrix\n", proyMat.toString());
  return proyMat;
}

/// Updates shader uniform with the MVP matrix
void updateMVPMatrix(ref Mat4f mvp, GLuint programId)
{
  // Get a handle for our "MVP" uniform
  // Only during the initialisation
  GLuint matrixID = glGetUniformLocation(programId, "MVP");
  // Send our transformation to the currently bound shader, in the "MVP" uniform
  // This is done in the main loop since each model will have a different MVP matrix (At least for the M part)
  // We need to transpose the matrix
  glUniformMatrix4fv(matrixID, 1, GL_TRUE, mvp.ptr);
}


void main()
{
  writeln("Keys 'a' and 'd' rotates the cube around Y axis. 'p' changes to perspective proyection, 'o' changes to ortho
      proyection. 'Esc' key exits the program");
  writeln("Press enter to continue");
  readln();

  scope(exit) {
    close();
  }
  if (!initSDL()) {
    return;
  }

  setGlFlags();
  if (!createAppWindow()) {
    return;
  }
  if (!initGl()) {
    return;
  }

  import std.math : PI_4;
  Qua_f cameraPos = Qua_f(Vec4f(0.3535534, 0.3535534, 0.1464466, 0.8535534)); // 3/4 view
  cameraPos.normalize();
  float modelAngle = 0;
  Vec3f modelPos = Vec3f.ZERO;
  bool ortho = true;

  auto programId = loadShaders();
  loadModel();
  updateScreen();

  // Polling for events
  bool quit = false;
  while(!quit) {
    SDL_PumpEvents();

    glUseProgram(programId);
    recalMatrixes(ortho, cameraPos, modelAngle, modelPos, programId);
    render();
    updateScreen();

    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      if (event.type == SDL_QUIT) {
        quit = true;
      }

      if (event.type == SDL_KEYDOWN) {
        import std.math : PI_4;
        const kevent = event.key;
        const keysym = kevent.keysym;
        switch (keysym.sym) {
          case SDL_Keycode.SDLK_ESCAPE:
          quit = true;
          break;

          case SDLK_o: // Toggle ortho mode
          ortho = true;
          break;

          case SDLK_p: // Toggle perspective mode
          ortho = false;
          break;

          case SDLK_a:
          modelAngle += PI_4 / 10.0f;
          break;

          case SDLK_d:
          modelAngle -= PI_4 / 10.0f;
          break;

          default:
        }
      }
    }
  }

}
