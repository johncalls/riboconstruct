#include <string>

#include <Patterncontainer.h>


/// Prototypes ////////////////////////////////////////////////////////////////
extern "C"
{
    char* evaluate(const char*, const char*, const char*, const char*);
    void  release(void);
}

char* p_pattern_str = NULL;


// hd_mode={n, s, e} (RNAf hairpin detection options)
// save_mode={c, d, r, i, s} (RNAf save mode options)

/// Implementation ////////////////////////////////////////////////////////////
char* evaluate(const char* p_hd_mode, const char* p_save_mode,
               const char* p_path, const char* p_param_file)
{
    std::string path(p_path);
    std::string param_file(p_param_file);
    Patterncontainer* p_pc = new Patterncontainer(path, param_file);
    if (!p_pc->readPatternFile())
        return NULL;

    for (unsigned int i = 0; i < p_pc->getNoOfPatterns(); ++i)
    {
        std::string hd_mode(p_hd_mode);
        std::string save_mode(p_save_mode);

        p_pc->executePattern(i, hd_mode);

        bool save_prog_call_info = false;
        for (unsigned int j = 0; j < save_mode.size(); ++j)
        {
            if (save_mode[j] == 'i')
                save_prog_call_info = true;
            else if (save_mode[j] == 's')
                p_pc->getPatternStr(i, &p_pattern_str);
        }
        ofstream file_info;
        if (save_prog_call_info)
        {
            std::string file;
            file.append(p_pc->getMainDir() + p_pc->getPatternName(i) + ".info");
            file_info.open(file.c_str());
            file_info << p_pc->listPattern(i);
            file_info.close();
        }

        if (save_mode.size() != 1 && save_mode[0] != 's')
            p_pc->savePattern(i, save_mode);
    }

    return p_pattern_str;
}


void release(void)
{
    if (p_pattern_str)
        delete p_pattern_str;
    p_pattern_str = NULL;
}
