/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GL2pdf_H_
#define model_GL2pdf_H_

#include <fstream>
#include <string>
#include <sstream>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#endif

#include <model/openGL/gl2ps.h>


namespace model {
	
    template <typename Renderer>
	class GL2pdf  {
		
        Renderer& renderer;
        
    public:
        
        static bool savePDF;
        
        
        GL2pdf(Renderer& r) :
        /* init list */ renderer(r)
        {/*! Constructor initializes renderer
          */
        }
        
        void saveAs(const std::string& filename)
        {/*! Saves to file
          */
            if(savePDF)
            {
                std::stringstream filenameStream;
                filenameStream << filename << ".pdf";
                std::string filenameWithExtension=filenameStream.str();
                std::cout<<"Saving file"<<filenameWithExtension<<std::endl;
                
                
                FILE *fp = fopen(filenameWithExtension.c_str(), "wb");
                GLint buffsize = 0,
                state = GL2PS_OVERFLOW;
                GLint viewport[4];
                
                glGetIntegerv(GL_VIEWPORT, viewport);
                
//                gl2psBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                gl2psBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                
                while( state == GL2PS_OVERFLOW ){
                    

                    
                    buffsize += 1024*1024;
                    gl2psBeginPage ( "MyTitle", "MySoftware", viewport,
//                                    GL2PS_EPS, GL2PS_BSP_SORT, GL2PS_SILENT |
                                    GL2PS_PDF, GL2PS_BSP_SORT, GL2PS_SILENT |
                                    GL2PS_SIMPLE_LINE_OFFSET |
//                                    GL2PS_SIMPLE_LINE_OFFSET | GL2PS_NO_BLENDING |
                                    GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT,
                                    GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
                                    fp, "MyFile" );
                    renderer.drawScene();
                    state = gl2psEndPage();
                }
                
                fclose(fp);
            }
            
        }
		
	};
	
    // Declare static data
    template <typename Renderer>
	bool GL2pdf<Renderer>::savePDF=true;
    
}
#endif
/*********************************************************************/
/*********************************************************************/





